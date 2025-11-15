// From standard library
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// core
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/foldertools.hpp>
#include <minimesh/core/util/numbers.hpp>
#include <minimesh/core/mohe/mesh_simplification.hpp>
#include <minimesh/core/mohe/laplacian_deformation.hpp>

// gui
#include <minimesh/viz/mesh_viewer.hpp>
#include <minimesh/viz/opengl_headers.hpp>




using namespace minimesh;

namespace globalvars
{
Mesh_viewer viewer;
mohecore::Mesh_connectivity mesh;
mohecore::Mesh_connectivity original_mesh;  // Store original mesh
//
int glut_main_window_id;
//
GLUI * glui;
//
int num_entities_to_simplify;
//
Eigen::Matrix3Xd displaced_vertex_positions;
Eigen::Matrix3Xd original_vertex_positions;

// Laplacian deformation variables
std::vector<int> fixed_handles;        // Indices of fixed vertices
std::vector<int> movable_handles;      // Indices of movable vertices
std::map<int, Eigen::Vector3d> handle_displacements;  // Displacement for each movable handle
int deformation_mode = 0;  // 0 = off, 1 = on (for GLUI checkbox)
int fix_vertex_mode = 0;   // 0 = off, 1 = on (for fixing vertices)
int selected_handle = -1;
int weight_choice = 0;  // 0 = cotangent, 1 = uniform

std::string subdivision_type;

}

// Forward declarations
namespace freeglutcallback {
	void update_handle_visualization();
}

// ======================================================
//              FREEGLUT CALL BACKS
// ======================================================
namespace freeglutcallback
{

void draw()
{
	globalvars::viewer.draw();
}


void window_reshaped(int w, int h)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.window_reshaped(w, h);

	if(should_redraw)
		glutPostRedisplay();
}


void keyboard_pressed(unsigned char c, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);

	if(should_redraw)
		glutPostRedisplay();
}


void keyboard_arrows_pressed(int c, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);

	if(should_redraw)
		glutPostRedisplay();
}

void mouse_pushed(int button, int state, int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.mouse_pushed(button, state, x, y);

	//
	// Handle vertex selection for fixing
	//
	{
		int clicked_on_vertex;
		bool did_user_click;
		globalvars::viewer.get_and_clear_vertex_selection(did_user_click, clicked_on_vertex);
		
		// If in fix vertex mode, add/remove vertex from fixed handles
		if(did_user_click && globalvars::fix_vertex_mode && clicked_on_vertex >= 0)
		{
			auto it = std::find(globalvars::fixed_handles.begin(), 
			                   globalvars::fixed_handles.end(), 
			                   clicked_on_vertex);
			if(it != globalvars::fixed_handles.end())
			{
				// Remove if already fixed
				globalvars::fixed_handles.erase(it);
				printf("\n[UNFIXED] Removed fixed vertex %d (total fixed: %zu)\n", 
				       clicked_on_vertex, globalvars::fixed_handles.size());
			}
			else
			{
				// Add to fixed handles
				globalvars::fixed_handles.push_back(clicked_on_vertex);
				
				// Get and print vertex position
				auto vert = globalvars::mesh.vertex_at(clicked_on_vertex);
				Eigen::Vector3d pos = vert.xyz();
				
				printf("\n=== FIXED VERTEX DEBUG INFO ===\n");
				printf("Vertex ID: %d\n", clicked_on_vertex);
				printf("Position: (%.6f, %.6f, %.6f)\n", pos.x(), pos.y(), pos.z());
				printf("Total fixed vertices: %zu\n", globalvars::fixed_handles.size());
				printf("Fixed vertex list: ");
				for(size_t i = 0; i < globalvars::fixed_handles.size(); ++i) {
					printf("%d", globalvars::fixed_handles[i]);
					if(i < globalvars::fixed_handles.size() - 1) printf(", ");
				}
				printf("\n==============================\n\n");
			}
			
			// Update visualization to show all fixed vertices
			freeglutcallback::update_handle_visualization();
			should_redraw = true;
		}
	}

	if(should_redraw)
		glutPostRedisplay();
}


void mouse_moved(int x, int y)
{
	bool should_redraw = false;
	should_redraw = should_redraw || globalvars::viewer.mouse_moved(x, y);

	{
		bool has_pull_performed;
		Eigen::Vector3f pull_amount;
		int pulled_vert;
		globalvars::viewer.get_and_clear_vertex_displacement(has_pull_performed, pull_amount, pulled_vert);

		if(has_pull_performed && globalvars::deformation_mode)
		{
			force_assert(pulled_vert != Mesh_viewer::invalid_index);

			try {
				mohecore::Laplacian_deformation deformer(globalvars::original_mesh);
				Eigen::MatrixXd positions(globalvars::original_mesh.n_active_vertices(), 3);
				Eigen::MatrixXi faces;
				
				int idx = 0;
				std::map<int, int> old_to_new;
				for(int vid = 0; vid < globalvars::original_mesh.n_total_vertices(); ++vid) {
					auto vert = globalvars::original_mesh.vertex_at(vid);
					if(vert.is_active()) {
						old_to_new[vid] = idx;
						positions.row(idx) = vert.xyz();
						idx++;
					}
				}
				
				std::vector<std::vector<int>> face_list;
				for(int fid = 0; fid < globalvars::original_mesh.n_total_faces(); ++fid) {
					auto face = globalvars::original_mesh.face_at(fid);
					if(face.is_active()) {
						auto he = face.half_edge();
						auto he_start = he;
						std::vector<int> fverts;
						do {
							fverts.push_back(old_to_new[he.origin().index()]);
							he = he.next();
						} while(!he.is_equal(he_start) && fverts.size() < 10);
						if(fverts.size() == 3) face_list.push_back(fverts);
					}
				}
				faces.resize(face_list.size(), 3);
				for(size_t i = 0; i < face_list.size(); ++i) {
					faces.row(i) << face_list[i][0], face_list[i][1], face_list[i][2];
				}
				
				std::vector<int> handles; 
				int n_constraints = globalvars::fixed_handles.size() + 1;
				Eigen::MatrixXd u_c(n_constraints, 3);
				
				for(size_t i = 0; i < globalvars::fixed_handles.size(); ++i) {
					int vid = globalvars::fixed_handles[i];
					handles.push_back(vid);
					u_c.row(i) = positions.row(vid);
				}
				
				handles.push_back(pulled_vert);
				auto current_vert = globalvars::mesh.vertex_at(pulled_vert);
				Eigen::RowVector3d current_deformed_pos = current_vert.xyz().transpose();
				Eigen::Vector3d pull_vec = pull_amount.cast<double>();
				Eigen::RowVector3d displacement = pull_vec.transpose();
				u_c.row(globalvars::fixed_handles.size()) = current_deformed_pos + displacement;
				
				auto p_final = deformer.solve_arap(positions, faces, handles, u_c, 3, 1e-4);
				
				idx = 0;
				for(int vid = 0; vid < globalvars::mesh.n_total_vertices(); ++vid) {
					auto vert = globalvars::mesh.vertex_at(vid);
					if(vert.is_active()) {
						vert.data().xyz = p_final.row(idx).transpose();
						globalvars::displaced_vertex_positions.col(idx) = p_final.row(idx).transpose();
						idx++;
					}
				}
				
				globalvars::viewer.get_mesh_buffer().set_vertex_positions(
					globalvars::displaced_vertex_positions.cast<float>());
				
				should_redraw = true;
			}
			catch(const std::exception& e) {
				printf("Deformation error: %s\n", e.what());
			}
		}
	}

	if(should_redraw)
		glutPostRedisplay();
}

void update_lowest_cost_edge_colors(){}


void subdivide_pressed(int)
{
	std::cout<<"not implemented yet"<<std::endl;
}


void simplify_pressed(int)
{

	std::cout<<"not implemented yet"<<std::endl;
}


void show_spheres_pressed(int)
{
  std::cout<<"not implemented yet"<<std::endl;
}


void update_handle_visualization()
{
	// Visualize handles with colored spheres
	int total_handles = globalvars::fixed_handles.size() + globalvars::movable_handles.size();
	
	if(total_handles == 0)
	{
		// Clear spheres - use a dummy vertex with transparent color
		Eigen::VectorXi sphere_indices(1);
		sphere_indices << 0;
		Eigen::Matrix4Xf sphere_colors(4, 1);
		sphere_colors << 0.0f, 0.0f, 0.0f, 0.0f;  // Fully transparent
		globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);
		return;
	}
	
	Eigen::VectorXi sphere_indices(total_handles);
	Eigen::Matrix4Xf sphere_colors(4, total_handles);
	
	int idx = 0;
	
	// Fixed handles - RED
	for(int v : globalvars::fixed_handles)
	{
		sphere_indices(idx) = v;
		sphere_colors.col(idx) << 1.0f, 0.0f, 0.0f, 1.0f;  // Red
		idx++;
	}
	
	// Movable handles - GREEN
	for(int v : globalvars::movable_handles)
	{
		sphere_indices(idx) = v;
		sphere_colors.col(idx) << 0.0f, 1.0f, 0.0f, 1.0f;  // Green
		idx++;
	}
	
	globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);
}

void fix_mode_callback(int)
{
	// "Fix" button - enters fix points mode (select vertices to fix)
	globalvars::fix_vertex_mode = 1;
	globalvars::deformation_mode = 0;
	globalvars::viewer.get_mouse_function() = Mesh_viewer::MOUSE_SELECT;
	printf("\n[MODE] Fix button clicked - FIX POINTS MODE\n");
	printf("  Click vertices to fix them (red spheres)\n");
	printf("  Click again to unfix\n");
	glutPostRedisplay();
}

void move_mode_callback(int)
{
	// "Move" button - enters drag points mode (deform by dragging)
	globalvars::deformation_mode = 1;
	globalvars::fix_vertex_mode = 0;
	globalvars::viewer.get_mouse_function() = Mesh_viewer::MOUSE_MOVE_VERTEX;
	printf("\n[MODE] Move button clicked - DRAG POINTS MODE\n");
	printf("  Drag vertices to deform the mesh\n");
	printf("  Fixed vertices (red) will stay in place\n");
	glutPostRedisplay();
}

void toggle_deformation_mode(int)
{
	globalvars::deformation_mode = !globalvars::deformation_mode;
	
	if(globalvars::deformation_mode)
	{
		printf("\n=== LAPLACIAN DEFORMATION MODE ENABLED ===\n");
		printf("Click vertices to add fixed handles (RED)\n");
		printf("Shift+Click vertices to add movable handles (GREEN)\n");
		printf("Click again to remove handles\n");
		printf("Use 'Move vertex' mode to drag movable handles\n");
		printf("==========================================\n\n");
		
		// Switch to select mode
		globalvars::viewer.get_mouse_function() = Mesh_viewer::MOUSE_SELECT;
	}
	else
	{
		printf("Deformation mode disabled\n");
	}
	
	freeglutcallback::update_handle_visualization();
	glutPostRedisplay();
}

void clear_handles_pressed(int)
{
	globalvars::fixed_handles.clear();
	globalvars::movable_handles.clear();
	globalvars::handle_displacements.clear();
	
	printf("Cleared all handles\n");
	
	freeglutcallback::update_handle_visualization();
	glutPostRedisplay();
}

void reset_mesh_pressed(int)
{
	// Restore original mesh
	globalvars::mesh.copy(globalvars::original_mesh);
	
	// Restore original positions
	globalvars::displaced_vertex_positions = globalvars::original_vertex_positions;
	
	// Clear handles and displacements
	globalvars::handle_displacements.clear();
	
	// Update viewer
	mohecore::Mesh_connectivity::Defragmentation_maps defrag;
	globalvars::mesh.compute_defragmention_maps(defrag);
	globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	
	// Update handle visualization
	freeglutcallback::update_handle_visualization();
	
	printf("Mesh reset to original state\n");
	glutPostRedisplay();
}

}



int main(int argc, char * argv[])
{
	// Remember current folder
	foldertools::pushd();
	// If no command line argument is specified, load a hardcoded mesh.
	// Useful when debugging with visual studio.
	// Change the hardcoded address to your needs.
	if(argc < 2)
	{
		throw std::runtime_error("Missing arguments meshfile or subdivision type.");
	}
	else // otherwise use the address specified in the command line
	{
    std::cout<<"hello minimesh!"<<std::endl;
		mohecore::Mesh_io(globalvars::mesh).read_auto(argv[1]);
		// globalvars::subdivision_type = argv[2];
	}

	// Initialize GLUT window
	glutInit(&argc, argv);
	glutInitWindowSize(800, 600);
	glutInitDisplayMode(GLUT_STENCIL | GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
	globalvars::glut_main_window_id = glutCreateWindow("Mesh Viewer");

	// Initialize GLUI window for buttons and ...
	globalvars::glui = GLUI_Master.create_glui("Controls");
	globalvars::glui->set_main_gfx_window(globalvars::glut_main_window_id);

	// Register callbacks
	glutDisplayFunc(freeglutcallback::draw);
	GLUI_Master.set_glutReshapeFunc(freeglutcallback::window_reshaped);
	GLUI_Master.set_glutKeyboardFunc(freeglutcallback::keyboard_pressed);
	GLUI_Master.set_glutSpecialFunc(freeglutcallback::keyboard_arrows_pressed);
	GLUI_Master.set_glutMouseFunc(freeglutcallback::mouse_pushed);
	glutMotionFunc(freeglutcallback::mouse_moved);
	GLUI_Master.set_glutIdleFunc(NULL);

	// Initialize the viewer (it needs the bounding box of the mesh)
	Eigen::AlignedBox3f bbox;
	for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
	{
		mohecore::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
		if(vertex.is_active())
		{
			bbox.extend(vertex.xyz().cast<float>());
		}
	}
	globalvars::viewer.initialize(bbox);

	// Load the mesh in the viewer
	{
		mohecore::Mesh_connectivity::Defragmentation_maps defrag;
		globalvars::mesh.compute_defragmention_maps(defrag);
		globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
	}
	
	// Removed edge collapse visualization for performance
	// freeglutcallback::update_lowest_cost_edge_colors();

	//
	// Add radio buttons to see which mesh components to view
	// Please view GLUI's user manual to learn more.
	//

	GLUI_Panel * panel_view = globalvars::glui->add_panel("View mesh components");
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show vertices", &globalvars::viewer.get_draw_vertices());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show edges", &globalvars::viewer.get_draw_edges());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show faces", &globalvars::viewer.get_draw_faces());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show axis", &globalvars::viewer.get_draw_axis());
	globalvars::glui->add_checkbox_to_panel(panel_view, "Show lighting", &globalvars::viewer.get_has_lighting());

	//
	// Add radio buttons to determine mouse left click functionality
	//
	GLUI_Panel * panel_mouse_func = globalvars::glui->add_panel("Mouse functionality");
	GLUI_RadioGroup * radio_group_mouse_func =   globalvars::glui->add_radiogroup_to_panel(panel_mouse_func, &globalvars::viewer.get_mouse_function());
	for(int i = 0; i < Mesh_viewer::MOUSE_INVALID; ++i)
	{
		if(i == Mesh_viewer::MOUSE_VIEW)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Pan and zoom");
		if(i == Mesh_viewer::MOUSE_SELECT)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Select vertex");
		if(i == Mesh_viewer::MOUSE_MOVE_VERTEX)
			globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Move vertex");
	}

	//
	// Add subdivide button
	//
	GLUI_Button* button_subdivide =  globalvars::glui->add_button("Subdivide Loop", -1, freeglutcallback::subdivide_pressed);
	button_subdivide->set_w(200);

	//
	// Add simplify button and a spinner to read how many entities to remove
	//
	globalvars::num_entities_to_simplify = 0;
	GLUI_Spinner* spinner_simplify = globalvars::glui->add_spinner("# of entities to simplify", GLUI_SPINNER_INT, &globalvars::num_entities_to_simplify);
	spinner_simplify->set_alignment(GLUI_ALIGN_CENTER);
	spinner_simplify->set_w(300);
	
	GLUI_Button* button_simplify = globalvars::glui->add_button("Simplify", -1, freeglutcallback::simplify_pressed);
	button_simplify->set_w(200);

	//
	// Add show spheres button to demo how to draw spheres on top of the vertices
	//
	globalvars::glui->add_button("Demo Showing Spheres", -1, freeglutcallback::show_spheres_pressed);


	//
	// Save the initial vertex positions and original mesh
	//
	globalvars::displaced_vertex_positions.resize(3, globalvars::mesh.n_active_vertices());
	globalvars::original_vertex_positions.resize(3, globalvars::mesh.n_active_vertices());
	for(int i = 0; i < globalvars::mesh.n_active_vertices(); ++i)
	{
		Eigen::Vector3d pos = globalvars::mesh.vertex_at(i).xyz();
		globalvars::displaced_vertex_positions.col(i) = pos;
		globalvars::original_vertex_positions.col(i) = pos;
	}
	
	// Save original mesh for reset
	globalvars::original_mesh.copy(globalvars::mesh);

	//
	// Add Laplacian Deformation controls
	//
	globalvars::glui->add_separator();
	
	GLUI_Panel* panel_laplacian = globalvars::glui->add_panel("Laplacian Deformation");
	
	// Fix button - enters fix points mode (select vertices to fix)
	GLUI_Button* btn_fix = globalvars::glui->add_button_to_panel(
		panel_laplacian, "Fix (click to fix)", -1,
		freeglutcallback::fix_mode_callback);
	btn_fix->set_w(250);
	
	// Move button - enters drag points mode (deform by dragging)
	GLUI_Button* btn_move = globalvars::glui->add_button_to_panel(
		panel_laplacian, "Move (drag to deform)", -1,
		freeglutcallback::move_mode_callback);
	btn_move->set_w(250);
	
	// Clear fixed vertices button
	GLUI_Button* button_clear_fixed = globalvars::glui->add_button_to_panel(
		panel_laplacian, 
		"Clear Fixed Vertices", 
		-1, 
		freeglutcallback::clear_handles_pressed);
	button_clear_fixed->set_w(250);
	
	globalvars::glui->add_separator_to_panel(panel_laplacian);
	
	// Weight type radio buttons
	GLUI_RadioGroup* radio_weight = globalvars::glui->add_radiogroup_to_panel(
		panel_laplacian, 
		&globalvars::weight_choice);
	globalvars::glui->add_radiobutton_to_group(radio_weight, "Cotangent weights");
	globalvars::glui->add_radiobutton_to_group(radio_weight, "Uniform weights");
	radio_weight->set_int_val(0);  // Default to cotangent
	
	globalvars::glui->add_separator_to_panel(panel_laplacian);
	
	// Reset mesh button
	GLUI_Button* button_reset = globalvars::glui->add_button_to_panel(
		panel_laplacian, 
		"Reset Mesh", 
		-1, 
		freeglutcallback::reset_mesh_pressed);
	button_reset->set_w(250);
	GLUI_StaticText* text5 = globalvars::glui->add_statictext("4. Switch to 'Move vertex' mode");
	GLUI_StaticText* text6 = globalvars::glui->add_statictext("5. Drag movable handles");
	GLUI_StaticText* text7 = globalvars::glui->add_statictext("6. Click 'Apply Deformation'");

	// Sync all glui variables
	globalvars::glui->sync_live();

	// Start main loop
	glutPostRedisplay(); // Draw everything again just for caution.
	glutMainLoop();

	// revert back to initial folder
	foldertools::popd();

	return 0;
}

