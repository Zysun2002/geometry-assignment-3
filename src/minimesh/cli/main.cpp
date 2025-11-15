//
// This is a bare executable that gets linked to the core library.
// You can use it to test your code, or start learning about the library.
//
// Here I have put an executable which reads a predefined mesh, flips a specific edge in that 
// mesh, and then writes it back to .vtk and .obj formats. Feel free to play with the example, change
// it, or move it to a completely different file.
//

#include <cstdio>
#include <iostream>

#define M_PI 3.14159265358979323846
#include <cmath>
#include <map>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;
fs::path out_path = "E:/Ziyu/workspace/course/geometryModeling/assignment_2/ZIyuSun_a2/mesh/output";
fs::path mesh_path = "E:/Ziyu/workspace/course/geometryModeling/assignment_2/ZIyuSun_a2/mesh";

#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/macros.hpp>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_simplification.hpp>
#include <minimesh/core/mohe/laplacian_deformation.hpp>

using namespace minimesh;

// ===================
// EXAMPLE UTILITY FUNCTIONS
// ===================

namespace 
{

//
// Create an example mesh file that we can read later.
//
void write_sphere_mesh()
{
    FILE *fl = fopen((out_path / "example_sphere.obj").string().c_str(), "w");
    if (!fl) return;

    // Vertices (icosahedron normalized)
    fprintf(fl, "v -0.525731 0.850651 0.000000\n");
    fprintf(fl, "v 0.525731 0.850651 0.000000\n");
    fprintf(fl, "v -0.525731 -0.850651 0.000000\n");
    fprintf(fl, "v 0.525731 -0.850651 0.000000\n");
    fprintf(fl, "v 0.000000 -0.525731 0.850651\n");
    fprintf(fl, "v 0.000000 0.525731 0.850651\n");
    fprintf(fl, "v 0.000000 -0.525731 -0.850651\n");
    fprintf(fl, "v 0.000000 0.525731 -0.850651\n");
    fprintf(fl, "v 0.850651 0.000000 -0.525731\n");
    fprintf(fl, "v 0.850651 0.000000 0.525731\n");
    fprintf(fl, "v -0.850651 0.000000 -0.525731\n");
    fprintf(fl, "v -0.850651 0.000000 0.525731\n\n");

    // Faces
    const int faces[20][3] = {
        {1,12,6}, {1,6,2}, {1,2,8}, {1,8,11}, {1,11,12},
        {2,6,10}, {6,12,5}, {12,11,3}, {11,8,7}, {8,2,9},
        {4,10,5}, {4,5,3}, {4,3,7}, {4,7,9}, {4,9,10},
        {5,10,6}, {3,5,12}, {7,3,11}, {9,7,8}, {10,9,2}
    };

    for (int i = 0; i < 20; ++i)
        fprintf(fl, "f %d %d %d\n", faces[i][0], faces[i][1], faces[i][2]);

    fclose(fl);
}

} // end of anonymus namespace

int main(int argc, char **argv)
{
  // Create mesh and Laplacian deformation object
  mohecore::Mesh_connectivity mesh;
  mohecore::Laplacian_deformation laplacian_def(mesh);

  // CLI functionality removed - use the GUI version instead
  std::cout << "CLI version is deprecated. Please use the GUI version (minimeshgui.exe)" << std::endl;

  return 0;
} // end of main()
