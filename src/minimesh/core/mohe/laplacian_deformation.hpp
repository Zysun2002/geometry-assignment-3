#pragma once

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
    Laplacian_deformation(Mesh_connectivity & mesh_in);

    Mesh_connectivity & grid() { return _m; }
    const Mesh_connectivity & grid() const { return _m; }

    Eigen::SparseMatrix<double> build_laplacian(
        const Eigen::MatrixXd& p,
        const Eigen::MatrixXi& faces);

    Eigen::MatrixXd compute_laplacian_coordinates(
        const Eigen::SparseMatrix<double>& L,
        const Eigen::MatrixXd& p);

    std::vector<Eigen::Matrix3d> compute_local_rotations(
        const Eigen::MatrixXd& p_original,
        const Eigen::MatrixXd& p_current,
        const Eigen::MatrixXi& faces);

    Eigen::MatrixXd solve_arap(
        const Eigen::MatrixXd& p_original,
        const Eigen::MatrixXi& faces,
        const std::vector<int>& handles,
        const Eigen::MatrixXd& u_c,
        int max_iterations = 10,
        double tolerance = 1e-6);

private:
    Mesh_connectivity & _m;

    double cotangent_weight(
        const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2);

    std::vector<int> get_free_indices(
        int n_vertices,
        const std::vector<int>& handles);
};

} // end of mohecore
} // end of minimesh
