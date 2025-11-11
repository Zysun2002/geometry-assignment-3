#pragma once

//
// laplacian_deformation.hpp
//
// Simplified Laplacian Mesh Deformation
// Implements linear solve-based deformation without rotation fitting
//

#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{

class Laplacian_deformation
{
public:
    // Constructor
    Laplacian_deformation(Mesh_connectivity & mesh_in);

    // Get the underlying mesh
    Mesh_connectivity & grid() { return _m; }
    const Mesh_connectivity & grid() const { return _m; }

    // ========================================================
    // Module 1: Build Laplacian matrix
    // ========================================================
    // Builds the Laplacian matrix L for the mesh
    // Input:
    //   p: Nx3 matrix of vertex positions (N = number of vertices)
    //   faces: Mx3 matrix of face indices (M = number of faces)
    //   weight_type: "cotangent" (default) or "uniform"
    // Output:
    //   L: NxN sparse Laplacian matrix
    Eigen::SparseMatrix<double> build_laplacian(
        const Eigen::MatrixXd& p,
        const Eigen::MatrixXi& faces,
        const std::string& weight_type = "cotangent");

    // ========================================================
    // Module 2: Compute Laplacian coordinates
    // ========================================================
    // Computes differential coordinates δ = L * p
    // Input:
    //   L: NxN sparse Laplacian matrix
    //   p: Nx3 matrix of vertex positions
    // Output:
    //   δ: Nx3 matrix of Laplacian coordinates
    Eigen::MatrixXd compute_laplacian_coordinates(
        const Eigen::SparseMatrix<double>& L,
        const Eigen::MatrixXd& p);

    // ========================================================
    // Module 3: Apply constraints
    // ========================================================
    // Splits the system into free and constrained parts
    // Input:
    //   L: NxN Laplacian matrix
    //   δ: Nx3 Laplacian coordinates
    //   handles: vector of constrained vertex indices
    //   u_c: |handles|x3 matrix of constraint positions
    // Output:
    //   L_ff: |free|x|free| submatrix for free vertices
    //   b_f: |free|x3 right-hand side (b_f = δ_f - L_fc * u_c)
    void apply_constraints(
        const Eigen::SparseMatrix<double>& L,
        const Eigen::MatrixXd& delta,
        const std::vector<int>& handles,
        const Eigen::MatrixXd& u_c,
        Eigen::SparseMatrix<double>& L_ff,
        Eigen::MatrixXd& b_f);

    // ========================================================
    // Module 4: Solve positions
    // ========================================================
    // Solves the linear system L_ff * p_f = b_f
    // Input:
    //   L_ff: |free|x|free| sparse matrix
    //   b_f: |free|x3 right-hand side
    // Output:
    //   p_f_solved: |free|x3 solution for free vertices
    Eigen::MatrixXd solve_positions(
        const Eigen::SparseMatrix<double>& L_ff,
        const Eigen::MatrixXd& b_f);

    // ========================================================
    // Module 5: Merge results
    // ========================================================
    // Combines free and constrained vertices into final result
    // Input:
    //   p: Nx3 original positions (for dimensions)
    //   free_idx: vector of free vertex indices
    //   cons_idx: vector of constrained vertex indices
    //   p_f_solved: |free|x3 solved free positions
    //   u_c: |cons|x3 constraint positions
    // Output:
    //   p_final: Nx3 final deformed positions
    Eigen::MatrixXd merge_results(
        const Eigen::MatrixXd& p,
        const std::vector<int>& free_idx,
        const std::vector<int>& cons_idx,
        const Eigen::MatrixXd& p_f_solved,
        const Eigen::MatrixXd& u_c);

    // ========================================================
    // Module 6: Run test
    // ========================================================
    // Complete test workflow
    // Loads woody-lo.obj, applies deformation, prints diagnostics
    void run_test();

    // ========================================================
    // ARAP: Compute local rotations
    // ========================================================
    // Computes optimal rotation matrix for each vertex based on neighbors
    // Input:
    //   p_original: Nx3 original vertex positions
    //   p_current: Nx3 current/deformed vertex positions
    //   faces: Mx3 face indices
    //   weight_type: "cotangent" or "uniform"
    // Output:
    //   rotations: vector of N 3x3 rotation matrices
    std::vector<Eigen::Matrix3d> compute_local_rotations(
        const Eigen::MatrixXd& p_original,
        const Eigen::MatrixXd& p_current,
        const Eigen::MatrixXi& faces,
        const std::string& weight_type = "cotangent");

    // ========================================================
    // ARAP: Iterative solver with rotation updates
    // ========================================================
    // Performs ARAP deformation with iterative rotation updates
    // Input:
    //   p_original: Nx3 original vertex positions
    //   faces: Mx3 face indices
    //   handles: vector of constrained vertex indices
    //   u_c: |handles|x3 constraint positions
    //   weight_type: "cotangent" or "uniform"
    //   max_iterations: maximum number of iterations (default 10)
    //   tolerance: convergence threshold (default 1e-6)
    // Output:
    //   p_final: Nx3 final deformed positions
    Eigen::MatrixXd solve_arap(
        const Eigen::MatrixXd& p_original,
        const Eigen::MatrixXi& faces,
        const std::vector<int>& handles,
        const Eigen::MatrixXd& u_c,
        const std::string& weight_type = "cotangent",
        int max_iterations = 10,
        double tolerance = 1e-6);

private:
    // Reference to the mesh
    Mesh_connectivity & _m;

    // Helper: Compute cotangent weight for an edge
    double cotangent_weight(
        const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2);

    // Helper: Extract vertex positions and faces from mesh
    void extract_mesh_data(
        Eigen::MatrixXd& positions,
        Eigen::MatrixXi& faces);

    // Helper: Update mesh with new positions
    void update_mesh_positions(const Eigen::MatrixXd& positions);

    // Helper: Get free vertex indices (complement of handles)
    std::vector<int> get_free_indices(
        int n_vertices,
        const std::vector<int>& handles);
};

} // end of mohecore
} // end of minimesh
