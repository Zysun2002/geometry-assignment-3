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


Laplacian_deformation::Laplacian_deformation(Mesh_connectivity & mesh_in)
: _m(mesh_in)
{
}

// heper: wriiten by ai
double Laplacian_deformation::cotangent_weight(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2)
{
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    double dot_product = edge1.dot(edge2);
    double cross_magnitude = edge1.cross(edge2).norm();
    
    if (cross_magnitude < 1e-10) {
        return 0.0;
    }
    
    return dot_product / cross_magnitude;
}

// helper: assist by ai
std::vector<int> Laplacian_deformation::get_free_indices(
    int n_vertices,
    const std::vector<int>& handles)
{
    std::set<int> constrained_vertices(handles.begin(), handles.end());

    std::vector<int> free_vertices;
    free_vertices.reserve(n_vertices - handles.size());
    
    for (int i = 0; i < n_vertices; ++i) {
        if (constrained_vertices.find(i) == constrained_vertices.end()) {
            free_vertices.push_back(i);
        }
    }
    
    return free_vertices;
}


// helper: assist by ai
Eigen::SparseMatrix<double> Laplacian_deformation::build_laplacian(
    const Eigen::MatrixXd& p,
    const Eigen::MatrixXi& faces)
{
    int num_vertices = p.rows();
    Eigen::SparseMatrix<double> L(num_vertices, num_vertices);
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(num_vertices * 8);
    
    std::vector<double> diagonal_sum(num_vertices, 0.0);
    
    // w_ij = (cot(alpha) + cot(beta)) 
    std::map<std::pair<int, int>, std::vector<int>> edge_to_triangles;
    
    // map edge to adjacent triangles
    for (int face_id = 0; face_id < faces.rows(); ++face_id) {
        for (int local_edge = 0; local_edge < 3; ++local_edge) {
            int v1 = faces(face_id, local_edge);
            int v2 = faces(face_id, (local_edge + 1) % 3);
            
            auto edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            edge_to_triangles[edge].push_back(face_id);
        }
    }
    
    // get cotangent weight
    for (const auto& [edge, adjacent_triangles] : edge_to_triangles) {
        int vertex_i = edge.first;
        int vertex_j = edge.second;
        
        double total_weight = 0.0;
        
        // sum cotangents from all triangles sharing this edge
        for (int triangle_id : adjacent_triangles) {
            // Find the vertex opposite to edge (i,j) in this triangle
            int opposite_vertex = -1;
            for (int corner = 0; corner < 3; ++corner) {
                int vertex = faces(triangle_id, corner);
                if (vertex != vertex_i && vertex != vertex_j) {
                    opposite_vertex = vertex;
                    break;
                }
            }
            
            if (opposite_vertex >= 0) {
                double cot_angle = cotangent_weight(
                    p.row(opposite_vertex),
                    p.row(vertex_i),
                    p.row(vertex_j)
                );
                total_weight += cot_angle;
            }
        }
        total_weight *= 0.5;
        
        // fill sparse matrix
        triplets.push_back(Eigen::Triplet<double>(vertex_i, vertex_j, -total_weight));
        triplets.push_back(Eigen::Triplet<double>(vertex_j, vertex_i, -total_weight));
        
        diagonal_sum[vertex_i] += total_weight;
        diagonal_sum[vertex_j] += total_weight;
    }
    
    // fill diagonal entries
    const double epsilon = 1e-8;
    for (int i = 0; i < num_vertices; ++i) {
        triplets.push_back(Eigen::Triplet<double>(i, i, diagonal_sum[i] + epsilon));
    }
    
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

// delta = L * p (prev)
Eigen::MatrixXd Laplacian_deformation::compute_laplacian_coordinates(
    const Eigen::SparseMatrix<double>& L,
    const Eigen::MatrixXd& p)
{
    return L * p;
}


// sovle rotation, assisted by ai
std::vector<Eigen::Matrix3d> Laplacian_deformation::compute_local_rotations(
    const Eigen::MatrixXd& p_original,
    const Eigen::MatrixXd& p_current,
    const Eigen::MatrixXi& faces)
{
    int num_vertices = p_original.rows();
    std::vector<Eigen::Matrix3d> rotations(num_vertices);
    
    // Store neighbors and cotangent weights for each vertex
    std::map<int, std::vector<std::pair<int, double>>> neighbors_and_weights;
    
    // prepare mapping data structure
    std::map<std::pair<int, int>, std::vector<int>> edge_to_triangles;
    
    for (int face_id = 0; face_id < faces.rows(); ++face_id) {
        for (int local_edge = 0; local_edge < 3; ++local_edge) {
            int v1 = faces(face_id, local_edge);
            int v2 = faces(face_id, (local_edge + 1) % 3);
            
            auto edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            edge_to_triangles[edge].push_back(face_id);
        }
    }
    
    // get cotangent weights again
    for (const auto& [edge, adjacent_triangles] : edge_to_triangles) {
        int vertex_i = edge.first;
        int vertex_j = edge.second;
        
        double total_weight = 0.0;
        for (int triangle_id : adjacent_triangles) {
            // buttferfly
            int opposite_vertex = -1;
            for (int corner = 0; corner < 3; ++corner) {
                int vertex = faces(triangle_id, corner);
                if (vertex != vertex_i && vertex != vertex_j) {
                    opposite_vertex = vertex;
                    break;
                }
            }
            
            if (opposite_vertex >= 0) {
                double cot_angle = cotangent_weight(
                    p_original.row(opposite_vertex),
                    p_original.row(vertex_i),
                    p_original.row(vertex_j)
                );
                total_weight += cot_angle;
            }
        }
        total_weight *= 0.5;
        
        neighbors_and_weights[vertex_i].push_back({vertex_j, total_weight});
        neighbors_and_weights[vertex_j].push_back({vertex_i, total_weight});
    }
    
    // get optimal rotation per vertex: assisted by ai
    for (int vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
        
        Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
        
        Eigen::RowVector3d pos_original = p_original.row(vertex_id);
        Eigen::RowVector3d pos_deformed = p_current.row(vertex_id);
        
        for (const auto& [neighbor_id, weight] : neighbors_and_weights[vertex_id]) {
            Eigen::Vector3d edge_original = (p_original.row(neighbor_id) - pos_original).transpose();
            Eigen::Vector3d edge_deformed = (p_current.row(neighbor_id) - pos_deformed).transpose();
            
            covariance += weight * edge_original * edge_deformed.transpose();
        }
        
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        
        Eigen::Matrix3d rotation = V * U.transpose();
        
        if (rotation.determinant() < 0) {
            Eigen::Matrix3d correction = Eigen::Matrix3d::Identity();
            correction(2, 2) = -1;
            rotation = V * correction * U.transpose();
        }
        
        rotations[vertex_id] = rotation;
    }
    
    return rotations;
}


// arap solver
Eigen::MatrixXd Laplacian_deformation::solve_arap(
    const Eigen::MatrixXd& p_original,
    const Eigen::MatrixXi& faces,
    const std::vector<int>& handles,
    const Eigen::MatrixXd& u_c,
    int max_iterations,
    double tolerance)
{
    int num_vertices = p_original.rows();
    
    Eigen::SparseMatrix<double> L = build_laplacian(p_original, faces);
    
    // very similart to buidling laplacian here
    std::map<int, std::vector<std::pair<int, double>>> neighbors_and_weights;
    std::map<std::pair<int, int>, std::vector<int>> edge_to_triangles;
    
    for (int face_id = 0; face_id < faces.rows(); ++face_id) {
        for (int local_edge = 0; local_edge < 3; ++local_edge) {
            int v1 = faces(face_id, local_edge);
            int v2 = faces(face_id, (local_edge + 1) % 3);
            auto edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            edge_to_triangles[edge].push_back(face_id);
        }
    }
    
    for (const auto& [edge, adjacent_triangles] : edge_to_triangles) {
        int vertex_i = edge.first;
        int vertex_j = edge.second;
        
        double total_weight = 0.0;
        for (int triangle_id : adjacent_triangles) {
            // butterfly 
            int opposite_vertex = -1;
            for (int corner = 0; corner < 3; ++corner) {
                int vertex = faces(triangle_id, corner);
                if (vertex != vertex_i && vertex != vertex_j) {
                    opposite_vertex = vertex;
                    break;
                }
            }
            if (opposite_vertex >= 0) {
                total_weight += cotangent_weight(
                    p_original.row(opposite_vertex),
                    p_original.row(vertex_i),
                    p_original.row(vertex_j)
                );
            }
        }
        total_weight *= 0.5;
        
        neighbors_and_weights[vertex_i].push_back({vertex_j, total_weight});
        neighbors_and_weights[vertex_j].push_back({vertex_i, total_weight});
    }
    
    std::vector<int> free_vertices = get_free_indices(num_vertices, handles);
    int num_free = free_vertices.size();
    int num_constrained = handles.size();
    
    //mappings: global vertex index -> local index in free/constrained arrays
    std::map<int, int> global_to_free, global_to_constrained;
    for (int i = 0; i < num_free; ++i) {
        global_to_free[free_vertices[i]] = i;
    }
    for (int i = 0; i < num_constrained; ++i) {
        global_to_constrained[handles[i]] = i;
    }
    
    // get only necessary blocks: L = [L_ff, L_fc; L_cf, L_cc]
    std::vector<Eigen::Triplet<double>> L_ff_triplets;
    std::vector<Eigen::Triplet<double>> L_fc_triplets;
    
    for (int col = 0; col < L.outerSize(); ++col) {
        // iterate through non-zero entries only
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it) {
            int row = it.row();
            int col_idx = it.col();
            double value = it.value();
            
            if (global_to_free.find(row) != global_to_free.end()) {
                int local_row = global_to_free[row];
                
                if (global_to_free.find(col_idx) != global_to_free.end()) {
                    // both are free
                    int local_col = global_to_free[col_idx];
                    L_ff_triplets.push_back(Eigen::Triplet<double>(local_row, local_col, value));
                } else if (global_to_constrained.find(col_idx) != global_to_constrained.end()) {
                    // one free, one constrained
                    int local_col = global_to_constrained[col_idx];
                    L_fc_triplets.push_back(Eigen::Triplet<double>(local_row, local_col, value));
                }
            }
        }
    }
    
    Eigen::SparseMatrix<double> L_ff(num_free, num_free);
    Eigen::SparseMatrix<double> L_fc(num_free, num_constrained);
    L_ff.setFromTriplets(L_ff_triplets.begin(), L_ff_triplets.end());
    L_fc.setFromTriplets(L_fc_triplets.begin(), L_fc_triplets.end());
    
    // precompute Cholesky factorization of L_ff
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky_solver;
    cholesky_solver.compute(L_ff);
    
    // solve first without rotations
    Eigen::MatrixXd delta_original = compute_laplacian_coordinates(L, p_original);
    Eigen::MatrixXd delta_free(num_free, 3);
    for (int i = 0; i < num_free; ++i) {
        delta_free.row(i) = delta_original.row(free_vertices[i]);
    }
    
    Eigen::MatrixXd target_delta = delta_free - L_fc * u_c;
    Eigen::MatrixXd p_current = p_original;
    
    Eigen::MatrixXd free_positions_solved(num_free, 3);
    for (int dim = 0; dim < 3; ++dim) {
        free_positions_solved.col(dim) = cholesky_solver.solve(target_delta.col(dim));
    }
    
    // assemble
    for (int i = 0; i < num_free; ++i) {
        p_current.row(free_vertices[i]) = free_positions_solved.row(i);
    }
    for (int i = 0; i < num_constrained; ++i) {
        p_current.row(handles[i]) = u_c.row(i);
    }
    
    // more complex iterative solver with rotations
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Eigen::MatrixXd p_previous = p_current;
        
        // get rotations
        std::vector<Eigen::Matrix3d> rotations = compute_local_rotations(p_original, p_current, faces);
        
        // with rotation fixed,
        Eigen::MatrixXd rotated_laplacian(num_vertices, 3);
        rotated_laplacian.setZero();
        
        for (int vertex_i = 0; vertex_i < num_vertices; ++vertex_i) {          
            for (const auto& [vertex_j, weight] : neighbors_and_weights[vertex_i]) {
                Eigen::Vector3d original_edge = (p_original.row(vertex_i) - p_original.row(vertex_j)).transpose();
                // eq 8: take average
                Eigen::Vector3d rotated_edge = 0.5 * (rotations[vertex_i] + rotations[vertex_j]) * original_edge;
                rotated_laplacian.row(vertex_i) += weight * rotated_edge.transpose();
            }
        }
        
        // prepare for solving
        Eigen::MatrixXd free_target(num_free, 3);
        for (int i = 0; i < num_free; ++i) {
            free_target.row(i) = rotated_laplacian.row(free_vertices[i]);
        }
        
        free_target -= L_fc * u_c;
        // solve it
        for (int dim = 0; dim < 3; ++dim) {
            free_positions_solved.col(dim) = cholesky_solver.solve(free_target.col(dim));
        }
        
        for (int i = 0; i < num_free; ++i) {
            p_current.row(free_vertices[i]) = free_positions_solved.row(i);
        }
        for (int i = 0; i < num_constrained; ++i) {
            p_current.row(handles[i]) = u_c.row(i);
        }
        
        // ai suggestion: early stopping
        double position_change = (p_current - p_previous).norm() / std::sqrt(num_vertices);
        if (position_change < tolerance) {
            break;
        }
    }
    
    return p_current;
}

} // end of mohecore
} // end of minimesh
