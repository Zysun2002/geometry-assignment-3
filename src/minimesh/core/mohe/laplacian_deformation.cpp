//
// laplacian_deformation.cpp
//
// Implementation of Simplified Laplacian Mesh Deformation
//

#include <minimesh/core/mohe/laplacian_deformation.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/util/assert.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <filesystem>

namespace fs = std::filesystem;

namespace minimesh
{
namespace mohecore
{

// Constructor
Laplacian_deformation::Laplacian_deformation(Mesh_connectivity & mesh_in)
: _m(mesh_in)
{
}

// ========================================================
// Helper: Compute cotangent weight
// ========================================================
double Laplacian_deformation::cotangent_weight(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2)
{
    // Compute cotangent of angle at v0
    // cot(θ) = cos(θ)/sin(θ) = (a·b)/(|a×b|)
    Eigen::Vector3d a = v1 - v0;
    Eigen::Vector3d b = v2 - v0;
    
    double dot = a.dot(b);
    double cross_norm = a.cross(b).norm();
    
    if (cross_norm < 1e-10) return 0.0; // Degenerate case
    
    return dot / cross_norm;
}

// ========================================================
// Helper: Extract mesh data
// ========================================================
void Laplacian_deformation::extract_mesh_data(
    Eigen::MatrixXd& positions,
    Eigen::MatrixXi& faces)
{
    int n_verts = _m.n_active_vertices();
    int n_faces = _m.n_active_faces();
    
    // Build mapping from mesh indices to continuous indices
    std::map<int, int> old_to_new;
    int new_idx = 0;
    
    positions.resize(n_verts, 3);
    
    // Extract vertex positions
    for (int vid = 0; vid < _m.n_total_vertices(); ++vid) {
        auto vert = _m.vertex_at(vid);
        if (vert.is_active()) {
            old_to_new[vid] = new_idx;
            positions.row(new_idx) = vert.xyz();
            new_idx++;
        }
    }
    
    // Extract face indices
    faces.resize(n_faces, 3);
    int face_idx = 0;
    
    for (int fid = 0; fid < _m.n_total_faces(); ++fid) {
        auto face = _m.face_at(fid);
        if (face.is_active()) {
            // Traverse face half-edges to get vertices
            auto he = face.half_edge();
            int vert_count = 0;
            std::vector<int> face_verts;
            
            auto he_start = he;
            do {
                face_verts.push_back(old_to_new[he.origin().index()]);
                he = he.next();
                vert_count++;
            } while (!he.is_equal(he_start) && vert_count < 10);
            
            // Only process triangular faces
            if (face_verts.size() == 3) {
                faces.row(face_idx) << face_verts[0], face_verts[1], face_verts[2];
                face_idx++;
            }
        }
    }
    
    // Resize if some faces were not triangular
    if (face_idx < n_faces) {
        faces.conservativeResize(face_idx, 3);
    }
}

// ========================================================
// Helper: Update mesh positions
// ========================================================
void Laplacian_deformation::update_mesh_positions(const Eigen::MatrixXd& positions)
{
    int new_idx = 0;
    for (int vid = 0; vid < _m.n_total_vertices(); ++vid) {
        auto vert = _m.vertex_at(vid);
        if (vert.is_active()) {
            vert.data().xyz = positions.row(new_idx);
            new_idx++;
        }
    }
}

// ========================================================
// Helper: Get free indices
// ========================================================
std::vector<int> Laplacian_deformation::get_free_indices(
    int n_vertices,
    const std::vector<int>& handles)
{
    // Create set of constrained indices for fast lookup
    std::map<int, bool> is_constrained;
    for (int h : handles) {
        is_constrained[h] = true;
    }
    
    // Collect free indices
    std::vector<int> free_idx;
    for (int i = 0; i < n_vertices; ++i) {
        if (is_constrained.find(i) == is_constrained.end()) {
            free_idx.push_back(i);
        }
    }
    
    return free_idx;
}

// ========================================================
// Module 1: Build Laplacian matrix
// ========================================================
Eigen::SparseMatrix<double> Laplacian_deformation::build_laplacian(
    const Eigen::MatrixXd& p,
    const Eigen::MatrixXi& faces,
    const std::string& weight_type)
{
    int n = p.rows();
    Eigen::SparseMatrix<double> L(n, n);
    
    // Use triplets for efficient construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * 8); // Estimate: each vertex has ~6-8 neighbors
    
    // Initialize diagonal sums
    std::vector<double> diag_sum(n, 0.0);
    
    if (weight_type == "cotangent") {
        // Cotangent weights
        // For each edge (i,j), weight w_ij = (cot α + cot β) / 2
        // where α and β are opposite angles
        
        // Build edge-to-triangles map
        std::map<std::pair<int, int>, std::vector<int>> edge_triangles;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                
                // Store edge with smaller index first
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                edge_triangles[edge].push_back(f);
            }
        }
        
        // Compute cotangent weights for each edge
        for (const auto& entry : edge_triangles) {
            int i = entry.first.first;
            int j = entry.first.second;
            const auto& tris = entry.second;
            
            double weight = 0.0;
            
            // Sum cotangents from all adjacent triangles
            for (int tri_idx : tris) {
                // Find the opposite vertex
                int k = -1;
                for (int v = 0; v < 3; ++v) {
                    int vid = faces(tri_idx, v);
                    if (vid != i && vid != j) {
                        k = vid;
                        break;
                    }
                }
                
                if (k >= 0) {
                    double cot = cotangent_weight(p.row(k), p.row(i), p.row(j));
                    weight += cot;
                }
            }
            
            // Average cotangent weight (divide by 2)
            weight *= 0.5;
            
            // Add to triplets (L is symmetric for undirected edges)
            triplets.push_back(Eigen::Triplet<double>(i, j, -weight));
            triplets.push_back(Eigen::Triplet<double>(j, i, -weight));
            
            diag_sum[i] += weight;
            diag_sum[j] += weight;
        }
    } else {
        // Uniform weights (weight = 1 for all edges)
        std::map<std::pair<int, int>, bool> edges;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                if (edges.find(edge) == edges.end()) {
                    edges[edge] = true;
                    
                    triplets.push_back(Eigen::Triplet<double>(i, j, -1.0));
                    triplets.push_back(Eigen::Triplet<double>(j, i, -1.0));
                    
                    diag_sum[i] += 1.0;
                    diag_sum[j] += 1.0;
                }
            }
        }
    }
    
    // Add diagonal entries (with regularization)
    const double regularization = 1e-8;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Eigen::Triplet<double>(i, i, diag_sum[i] + regularization));
    }
    
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

// ========================================================
// Module 2: Compute Laplacian coordinates
// ========================================================
Eigen::MatrixXd Laplacian_deformation::compute_laplacian_coordinates(
    const Eigen::SparseMatrix<double>& L,
    const Eigen::MatrixXd& p)
{
    // δ = L * p
    // Compute for each coordinate (x, y, z) independently
    Eigen::MatrixXd delta = L * p;
    return delta;
}

// ========================================================
// Module 3: Apply constraints
// ========================================================
void Laplacian_deformation::apply_constraints(
    const Eigen::SparseMatrix<double>& L,
    const Eigen::MatrixXd& delta,
    const std::vector<int>& handles,
    const Eigen::MatrixXd& u_c,
    Eigen::SparseMatrix<double>& L_ff,
    Eigen::MatrixXd& b_f)
{
    int n = L.rows();
    
    // Get free indices
    std::vector<int> free_idx = get_free_indices(n, handles);
    int n_free = free_idx.size();
    int n_cons = handles.size();
    
    // Build index mappings
    std::map<int, int> free_map; // old index -> new free index
    std::map<int, int> cons_map; // old index -> new constrained index
    
    for (int i = 0; i < n_free; ++i) {
        free_map[free_idx[i]] = i;
    }
    for (int i = 0; i < n_cons; ++i) {
        cons_map[handles[i]] = i;
    }
    
    // Extract L_ff (free-to-free) and L_fc (free-to-constrained)
    std::vector<Eigen::Triplet<double>> triplets_ff;
    Eigen::SparseMatrix<double> L_fc(n_free, n_cons);
    std::vector<Eigen::Triplet<double>> triplets_fc;
    
    // Iterate over non-zero entries of L
    for (int k = 0; k < L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double val = it.value();
            
            // Check if row is free
            if (free_map.find(i) != free_map.end()) {
                int i_new = free_map[i];
                
                // Check if column is free
                if (free_map.find(j) != free_map.end()) {
                    int j_new = free_map[j];
                    triplets_ff.push_back(Eigen::Triplet<double>(i_new, j_new, val));
                }
                // Check if column is constrained
                else if (cons_map.find(j) != cons_map.end()) {
                    int j_new = cons_map[j];
                    triplets_fc.push_back(Eigen::Triplet<double>(i_new, j_new, val));
                }
            }
        }
    }
    
    // Build L_ff
    L_ff.resize(n_free, n_free);
    L_ff.setFromTriplets(triplets_ff.begin(), triplets_ff.end());
    
    // Build L_fc
    L_fc.setFromTriplets(triplets_fc.begin(), triplets_fc.end());
    
    // Extract δ_f (Laplacian coordinates for free vertices)
    Eigen::MatrixXd delta_f(n_free, 3);
    for (int i = 0; i < n_free; ++i) {
        delta_f.row(i) = delta.row(free_idx[i]);
    }
    
    // Compute right-hand side: b_f = δ_f - L_fc * u_c
    b_f = delta_f - L_fc * u_c;
}

// ========================================================
// Module 4: Solve positions
// ========================================================
Eigen::MatrixXd Laplacian_deformation::solve_positions(
    const Eigen::SparseMatrix<double>& L_ff,
    const Eigen::MatrixXd& b_f)
{
    // Solve L_ff * p_f = b_f using SimplicialLDLT (for symmetric positive definite systems)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    
    // Call compute() only once per matrix (this is the expensive operation)
    solver.compute(L_ff);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERROR: SimplicialLDLT decomposition failed!" << std::endl;
        return Eigen::MatrixXd::Zero(b_f.rows(), b_f.cols());
    }
    
    // Solve for each coordinate independently
    // solve() can be called multiple times efficiently
    Eigen::MatrixXd p_f_solved(b_f.rows(), 3);
    
    for (int dim = 0; dim < 3; ++dim) {
        p_f_solved.col(dim) = solver.solve(b_f.col(dim));
        
        if (solver.info() != Eigen::Success) {
            std::cerr << "ERROR: Sparse solve failed for dimension " << dim << std::endl;
        }
    }
    
    return p_f_solved;
}

// ========================================================
// Module 5: Merge results
// ========================================================
Eigen::MatrixXd Laplacian_deformation::merge_results(
    const Eigen::MatrixXd& p,
    const std::vector<int>& free_idx,
    const std::vector<int>& cons_idx,
    const Eigen::MatrixXd& p_f_solved,
    const Eigen::MatrixXd& u_c)
{
    int n = p.rows();
    Eigen::MatrixXd p_final = p; // Copy original positions
    
    // Fill in free vertices
    for (size_t i = 0; i < free_idx.size(); ++i) {
        p_final.row(free_idx[i]) = p_f_solved.row(i);
    }
    
    // Fill in constrained vertices
    for (size_t i = 0; i < cons_idx.size(); ++i) {
        p_final.row(cons_idx[i]) = u_c.row(i);
    }
    
    return p_final;
}

// ========================================================
// ARAP: Compute local rotations
// ========================================================
std::vector<Eigen::Matrix3d> Laplacian_deformation::compute_local_rotations(
    const Eigen::MatrixXd& p_original,
    const Eigen::MatrixXd& p_current,
    const Eigen::MatrixXi& faces,
    const std::string& weight_type)
{
    int n = p_original.rows();
    std::vector<Eigen::Matrix3d> rotations(n);
    
    // Build edge weights and neighbor information
    std::map<int, std::vector<std::pair<int, double>>> vertex_neighbors;
    
    if (weight_type == "cotangent") {
        // Build edge-to-triangles map for cotangent weights
        std::map<std::pair<int, int>, std::vector<int>> edge_triangles;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                edge_triangles[edge].push_back(f);
            }
        }
        
        // Compute cotangent weights
        for (const auto& entry : edge_triangles) {
            int i = entry.first.first;
            int j = entry.first.second;
            const auto& tris = entry.second;
            
            double weight = 0.0;
            for (int tri_idx : tris) {
                int k = -1;
                for (int v = 0; v < 3; ++v) {
                    int vid = faces(tri_idx, v);
                    if (vid != i && vid != j) {
                        k = vid;
                        break;
                    }
                }
                
                if (k >= 0) {
                    double cot = cotangent_weight(p_original.row(k), p_original.row(i), p_original.row(j));
                    weight += cot;
                }
            }
            weight *= 0.5;
            
            vertex_neighbors[i].push_back({j, weight});
            vertex_neighbors[j].push_back({i, weight});
        }
    } else {
        // Uniform weights
        std::map<std::pair<int, int>, bool> edges;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                if (edges.find(edge) == edges.end()) {
                    edges[edge] = true;
                    vertex_neighbors[i].push_back({j, 1.0});
                    vertex_neighbors[j].push_back({i, 1.0});
                }
            }
        }
    }
    
    // Compute rotation for each vertex
    for (int i = 0; i < n; ++i) {
        if (vertex_neighbors.find(i) == vertex_neighbors.end() || vertex_neighbors[i].empty()) {
            rotations[i] = Eigen::Matrix3d::Identity();
            continue;
        }
        
        // Build covariance matrix using Kabsch algorithm formula
        // S = Σ_j w_ij * e_orig * e_curr^T (original * current^T)
        // This finds rotation R that best maps original edges to current deformed edges
        Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
        
        Eigen::RowVector3d pi_orig = p_original.row(i);
        Eigen::RowVector3d pi_curr = p_current.row(i);
        
        for (const auto& neighbor : vertex_neighbors[i]) {
            int j = neighbor.first;
            double w_ij = neighbor.second;
            
            Eigen::Vector3d e_orig = (p_original.row(j) - pi_orig).transpose();  // Column vector
            Eigen::Vector3d e_curr = (p_current.row(j) - pi_curr).transpose();   // Column vector
            
            // Outer product: (3x1) * (1x3) = (3x3)
            // Kabsch: S = Σ w_ij * e_orig * e_curr^T
            S += w_ij * e_orig * e_curr.transpose();
        }
        
        // Compute SVD: S = U * Σ * V^T
        // For Kabsch algorithm with S = Σ w_ij * e_orig * e_curr^T
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        Eigen::Vector3d singular_values = svd.singularValues();
        
        // Optimal rotation: R = V * U^T (Kabsch formula)
        // This rotation maps original edges to current edges
        Eigen::Matrix3d R = V * U.transpose();
        
        double det_R = R.determinant();
        
        // Ensure proper rotation (det(R) = 1, not -1 which would be a reflection)
        // For 2D meshes (planar, z=0), the third singular value is zero
        // We need to handle the reflection correction carefully
        if (det_R < 0) {
            // Create a diagonal matrix to fix the reflection
            // Flip the component corresponding to the smallest singular value
            Eigen::Matrix3d correction = Eigen::Matrix3d::Identity();
            correction(2, 2) = -1;  // Flip the third axis (smallest singular value for 2D meshes)
            R = V * correction * U.transpose();
        }
        
        // Debug output for first few vertices
        if (i < 3) {
            std::cout << "Vertex " << i << " rotation:" << std::endl;
            std::cout << "  Covariance S norm: " << S.norm() << std::endl;
            std::cout << "  Singular values: " << singular_values.transpose() << std::endl;
            std::cout << "  Det before fix: " << det_R << ", after: " << R.determinant() << std::endl;
            std::cout << "  R matrix:" << std::endl << R << std::endl;
            std::cout << "  R(2,2) component: " << R(2,2) << std::endl;
        }
        
        rotations[i] = R;
    }
    
    return rotations;
}

// ========================================================
// ARAP: Iterative solver with rotation updates
// ========================================================
Eigen::MatrixXd Laplacian_deformation::solve_arap(
    const Eigen::MatrixXd& p_original,
    const Eigen::MatrixXi& faces,
    const std::vector<int>& handles,
    const Eigen::MatrixXd& u_c,
    const std::string& weight_type,
    int max_iterations,
    double tolerance)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "ARAP SOLVER CALLED" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Vertices: " << p_original.rows() << std::endl;
    std::cout << "Handles: " << handles.size() << std::endl;
    std::cout << "Max iterations: " << max_iterations << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    int n = p_original.rows();
    
    // Build Laplacian matrix (only once)
    Eigen::SparseMatrix<double> L = build_laplacian(p_original, faces, weight_type);
    
    // Build edge weights and neighbor information for RHS computation
    std::map<int, std::vector<std::pair<int, double>>> vertex_neighbors;
    
    if (weight_type == "cotangent") {
        std::map<std::pair<int, int>, std::vector<int>> edge_triangles;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                edge_triangles[edge].push_back(f);
            }
        }
        
        for (const auto& entry : edge_triangles) {
            int i = entry.first.first;
            int j = entry.first.second;
            const auto& tris = entry.second;
            
            double weight = 0.0;
            for (int tri_idx : tris) {
                int k = -1;
                for (int v = 0; v < 3; ++v) {
                    int vid = faces(tri_idx, v);
                    if (vid != i && vid != j) {
                        k = vid;
                        break;
                    }
                }
                if (k >= 0) {
                    weight += cotangent_weight(p_original.row(k), p_original.row(i), p_original.row(j));
                }
            }
            weight *= 0.5;
            
            vertex_neighbors[i].push_back({j, weight});
            vertex_neighbors[j].push_back({i, weight});
        }
    } else {
        std::map<std::pair<int, int>, bool> edges;
        
        for (int f = 0; f < faces.rows(); ++f) {
            for (int k = 0; k < 3; ++k) {
                int i = faces(f, k);
                int j = faces(f, (k + 1) % 3);
                auto edge = std::make_pair(std::min(i, j), std::max(i, j));
                if (edges.find(edge) == edges.end()) {
                    edges[edge] = true;
                    vertex_neighbors[i].push_back({j, 1.0});
                    vertex_neighbors[j].push_back({i, 1.0});
                }
            }
        }
    }
    
    // Get free vertex indices
    std::vector<int> free_idx = get_free_indices(n, handles);
    int n_free = free_idx.size();
    int n_cons = handles.size();
    
    // Build index mappings
    std::map<int, int> free_map, cons_map;
    for (int i = 0; i < n_free; ++i) {
        free_map[free_idx[i]] = i;
    }
    for (int i = 0; i < n_cons; ++i) {
        cons_map[handles[i]] = i;
    }
    
    // Extract L_ff and L_fc (only once)
    std::vector<Eigen::Triplet<double>> triplets_ff;
    Eigen::SparseMatrix<double> L_fc(n_free, n_cons);
    std::vector<Eigen::Triplet<double>> triplets_fc;
    
    for (int k = 0; k < L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double val = it.value();
            
            if (free_map.find(i) != free_map.end()) {
                int i_new = free_map[i];
                if (free_map.find(j) != free_map.end()) {
                    int j_new = free_map[j];
                    triplets_ff.push_back(Eigen::Triplet<double>(i_new, j_new, val));
                } else if (cons_map.find(j) != cons_map.end()) {
                    int j_new = cons_map[j];
                    triplets_fc.push_back(Eigen::Triplet<double>(i_new, j_new, val));
                }
            }
        }
    }
    
    Eigen::SparseMatrix<double> L_ff(n_free, n_free);
    L_ff.setFromTriplets(triplets_ff.begin(), triplets_ff.end());
    L_fc.setFromTriplets(triplets_fc.begin(), triplets_fc.end());
    
    // Pre-factor L_ff (only once - expensive operation)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_ff);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERROR: SimplicialLDLT decomposition failed!" << std::endl;
        return p_original;
    }
    
    // Initialize with single linear solve (no rotations)
    Eigen::MatrixXd delta = compute_laplacian_coordinates(L, p_original);
    Eigen::MatrixXd delta_f(n_free, 3);
    for (int i = 0; i < n_free; ++i) {
        delta_f.row(i) = delta.row(free_idx[i]);
    }
    
    Eigen::MatrixXd b_f = delta_f - L_fc * u_c;
    Eigen::MatrixXd p_current = p_original;
    
    // Solve initial positions
    Eigen::MatrixXd p_f_solved(n_free, 3);
    for (int dim = 0; dim < 3; ++dim) {
        p_f_solved.col(dim) = solver.solve(b_f.col(dim));
    }
    
    // Merge initial solution
    for (int i = 0; i < n_free; ++i) {
        p_current.row(free_idx[i]) = p_f_solved.row(i);
    }
    for (int i = 0; i < n_cons; ++i) {
        p_current.row(handles[i]) = u_c.row(i);
    }
    
    // Check initial solution sanity
    double init_change = (p_current - p_original).norm() / std::sqrt(n);
    std::cout << "Initial linear solve change: " << init_change << std::endl;
    
    // ARAP iterations
    std::cout << "Starting ARAP iterations..." << std::endl;
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::MatrixXd p_previous = p_current;
        
        std::cout << "\n--- Iteration " << (iter + 1) << " ---" << std::endl;
        
        // Step 1: Update rotations
        std::vector<Eigen::Matrix3d> rotations = compute_local_rotations(p_original, p_current, faces, weight_type);
        
        // Step 2: Update right-hand side with rotations
        // For each vertex i: b_i = Σ_j w_ij * (R_i + R_j)/2 * (p_j^0 - p_i^0)
        // This represents the rotated differential coordinates
        Eigen::MatrixXd b_arap(n, 3);
        b_arap.setZero();
        
        for (int i = 0; i < n; ++i) {
            if (vertex_neighbors.find(i) == vertex_neighbors.end()) continue;
            
            for (const auto& neighbor : vertex_neighbors[i]) {
                int j = neighbor.first;
                double w_ij = neighbor.second;
                
                // Original edge from i to j
                Eigen::Vector3d e_ij = (p_original.row(i) - p_original.row(j)).transpose();

                
                // Apply average rotation: (R_i + R_j)/2 * e_ij
                Eigen::Vector3d rotated_edge = 0.5 * (rotations[i] + rotations[j]) * e_ij;
                
                // Accumulate weighted rotated edge
                b_arap.row(i) += w_ij * rotated_edge.transpose();
            }
        }
        
        // Extract b_f for free vertices
        Eigen::MatrixXd b_f_arap(n_free, 3);
        for (int i = 0; i < n_free; ++i) {
            b_f_arap.row(i) = b_arap.row(free_idx[i]);
        }
        
        // Subtract constrained contribution
        b_f_arap -= L_fc * u_c;
        
        // Step 3: Solve for new positions
        for (int dim = 0; dim < 3; ++dim) {
            p_f_solved.col(dim) = solver.solve(b_f_arap.col(dim));
        }
        
        // Merge solution
        for (int i = 0; i < n_free; ++i) {
            p_current.row(free_idx[i]) = p_f_solved.row(i);
        }
        for (int i = 0; i < n_cons; ++i) {
            p_current.row(handles[i]) = u_c.row(i);
        }
        
        // Check convergence
        double change = (p_current - p_previous).norm() / std::sqrt(n);
        
        // Debug: track specific vertex positions
        if (iter < 3) {
            std::cout << "  First vertex position: " << p_current.row(0) << std::endl;
            std::cout << "  Position change: " << (p_current.row(0) - p_previous.row(0)) << std::endl;
        }
        
        std::cout << "  Iteration " << (iter + 1) << ": change = " << change << std::endl;
        
        if (change < tolerance) {
            std::cout << "ARAP converged after " << (iter + 1) << " iterations" << std::endl;
            break;
        }
        
        if (iter == max_iterations - 1) {
            std::cout << "ARAP reached max iterations (" << max_iterations << ")" << std::endl;
        }
    }
    
    return p_current;
}

// ========================================================
// Module 6: Run test
// ========================================================
void Laplacian_deformation::run_test()
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Laplacian Mesh Deformation Test" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Step 1: Load mesh
    std::cout << "Step 1: Loading woody-lo.obj..." << std::endl;
    
    // Try multiple possible paths
    fs::path mesh_path;
    std::vector<fs::path> possible_paths = {
        fs::current_path() / "meshes" / "woody-lo.obj",
        fs::current_path().parent_path() / "meshes" / "woody-lo.obj",
        fs::current_path().parent_path().parent_path() / "meshes" / "woody-lo.obj",
        "meshes/woody-lo.obj",
        "../meshes/woody-lo.obj",
        "../../meshes/woody-lo.obj"
    };
    
    bool found = false;
    for (const auto& path : possible_paths) {
        if (fs::exists(path)) {
            mesh_path = path;
            found = true;
            break;
        }
    }
    
    if (!found) {
        std::cerr << "ERROR: Could not find woody-lo.obj in any expected location." << std::endl;
        std::cerr << "Current directory: " << fs::current_path() << std::endl;
        return;
    }
    
    std::cout << "  Found mesh at: " << mesh_path << std::endl;
    
    Mesh_io io(_m);
    io.read_auto(mesh_path.string());
    
    std::cout << "  Vertices: " << _m.n_active_vertices() << std::endl;
    std::cout << "  Faces: " << _m.n_active_faces() << std::endl;
    
    // Step 2: Extract mesh data
    std::cout << "\nStep 2: Extracting mesh data..." << std::endl;
    Eigen::MatrixXd p;
    Eigen::MatrixXi faces;
    extract_mesh_data(p, faces);
    
    std::cout << "  Position matrix: " << p.rows() << " x " << p.cols() << std::endl;
    std::cout << "  Face matrix: " << faces.rows() << " x " << faces.cols() << std::endl;
    
    // Step 3: Build Laplacian
    std::cout << "\nStep 3: Building Laplacian matrix with cotangent weights..." << std::endl;
    Eigen::SparseMatrix<double> L = build_laplacian(p, faces, "cotangent");
    std::cout << "  Laplacian size: " << L.rows() << " x " << L.cols() << std::endl;
    std::cout << "  Non-zeros: " << L.nonZeros() << std::endl;
    
    // Step 4: Compute Laplacian coordinates
    std::cout << "\nStep 4: Computing Laplacian coordinates..." << std::endl;
    Eigen::MatrixXd delta = compute_laplacian_coordinates(L, p);
    std::cout << "  Delta shape: " << delta.rows() << " x " << delta.cols() << std::endl;
    
    // Step 5: Define constraints (handle vertices)
    std::cout << "\nStep 5: Defining constraints..." << std::endl;
    
    // Find vertices to constrain (e.g., bottom vertices and a top vertex to pull)
    std::vector<int> handles;
    Eigen::VectorXd y_coords = p.col(1); // Y coordinates
    double y_min = y_coords.minCoeff();
    double y_max = y_coords.maxCoeff();
    double y_range = y_max - y_min;
    
    // Constrain bottom vertices (y < y_min + 0.1 * range)
    std::vector<int> bottom_handles;
    for (int i = 0; i < p.rows(); ++i) {
        if (p(i, 1) < y_min + 0.1 * y_range) {
            bottom_handles.push_back(i);
        }
    }
    
    // Find a top vertex to pull (y > y_max - 0.1 * range)
    std::vector<int> top_handles;
    for (int i = 0; i < p.rows(); ++i) {
        if (p(i, 1) > y_max - 0.1 * y_range) {
            top_handles.push_back(i);
            if (top_handles.size() >= 3) break; // Just a few top vertices
        }
    }
    
    handles.insert(handles.end(), bottom_handles.begin(), bottom_handles.end());
    handles.insert(handles.end(), top_handles.begin(), top_handles.end());
    
    std::cout << "  Bottom handles: " << bottom_handles.size() << std::endl;
    std::cout << "  Top handles: " << top_handles.size() << std::endl;
    std::cout << "  Total handles: " << handles.size() << std::endl;
    
    // Define constraint positions
    Eigen::MatrixXd u_c(handles.size(), 3);
    
    // Bottom stays fixed
    for (size_t i = 0; i < bottom_handles.size(); ++i) {
        u_c.row(i) = p.row(bottom_handles[i]);
    }
    
    // Top gets pulled upward and sideways
    for (size_t i = 0; i < top_handles.size(); ++i) {
        u_c.row(bottom_handles.size() + i) = p.row(top_handles[i]);
        u_c(bottom_handles.size() + i, 0) += 0.5; // Pull in X
        u_c(bottom_handles.size() + i, 1) += 1.0; // Pull in Y (up)
    }
    
    // Step 6: Apply constraints
    std::cout << "\nStep 6: Applying constraints..." << std::endl;
    Eigen::SparseMatrix<double> L_ff;
    Eigen::MatrixXd b_f;
    apply_constraints(L, delta, handles, u_c, L_ff, b_f);
    
    std::vector<int> free_idx = get_free_indices(p.rows(), handles);
    std::cout << "  Free vertices: " << free_idx.size() << std::endl;
    std::cout << "  L_ff size: " << L_ff.rows() << " x " << L_ff.cols() << std::endl;
    std::cout << "  b_f shape: " << b_f.rows() << " x " << b_f.cols() << std::endl;
    
    // Step 7: Solve for free positions
    std::cout << "\nStep 7: Solving linear system..." << std::endl;
    Eigen::MatrixXd p_f_solved = solve_positions(L_ff, b_f);
    std::cout << "  Solution shape: " << p_f_solved.rows() << " x " << p_f_solved.cols() << std::endl;
    
    // Step 8: Merge results
    std::cout << "\nStep 8: Merging results..." << std::endl;
    Eigen::MatrixXd p_final = merge_results(p, free_idx, handles, p_f_solved, u_c);
    std::cout << "  Final positions: " << p_final.rows() << " x " << p_final.cols() << std::endl;
    
    // Step 9: Compute diagnostics
    std::cout << "\n========================================" << std::endl;
    std::cout << "Diagnostics" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Constraint error: ||p_final[handles] - u_c||
    double constraint_error = 0.0;
    for (size_t i = 0; i < handles.size(); ++i) {
        constraint_error += (p_final.row(handles[i]) - u_c.row(i)).squaredNorm();
    }
    constraint_error = std::sqrt(constraint_error / handles.size());
    std::cout << "Constraint error (RMSE): " << constraint_error << std::endl;
    
    // Residual norm: ||L_ff * p_f_solved - b_f||
    Eigen::MatrixXd residual = L_ff * p_f_solved - b_f;
    double residual_norm = residual.norm() / std::sqrt(residual.size());
    std::cout << "Residual norm (RMS): " << residual_norm << std::endl;
    
    // Print some final positions
    std::cout << "\nFinal positions (first 5 vertices):" << std::endl;
    for (int i = 0; i < std::min(5, (int)p_final.rows()); ++i) {
        std::cout << "  v" << i << ": ["
                  << p_final(i, 0) << ", "
                  << p_final(i, 1) << ", "
                  << p_final(i, 2) << "]" << std::endl;
    }
    
    // Step 10: Update mesh and save
    std::cout << "\nStep 9: Updating mesh and saving..." << std::endl;
    update_mesh_positions(p_final);
    
    // Try to save in multiple locations
    std::vector<fs::path> out_paths = {
        fs::current_path() / "deformed_woody.obj",
        fs::current_path() / "mesh" / "deformed_woody.obj",
        "deformed_woody.obj"
    };
    
    fs::path out_path = out_paths[0];
    try {
        // Create directory if needed
        if (out_path.has_parent_path()) {
            fs::create_directories(out_path.parent_path());
        }
        io.write_obj(out_path.string());
        std::cout << "  Saved to: " << out_path << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Warning: Could not save to " << out_path << std::endl;
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test completed successfully!" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

} // end of mohecore
} // end of minimesh
