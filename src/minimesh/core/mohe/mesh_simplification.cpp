#include <minimesh/core/mohe/mesh_simplification.hpp>
#include <minimesh/core/util/assert.hpp>

#include <iostream>
#include <map>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

#define M_PI 3.14159265358979323846


namespace minimesh
{
namespace mohecore
{

Quadric::Quadric() : Q(Eigen::Matrix4d::Zero()) {}

Quadric::Quadric(const Eigen::Vector4d& plane) : Q(Eigen::Matrix4d::Zero()) 
{
    // Q = plane * plane^T
    Q = plane * plane.transpose();
}

Quadric Quadric::operator+(const Quadric& other) const 
{
    Quadric result;
    result.Q = Q + other.Q;
    return result;
}

double Quadric::computeError(const Eigen::Vector3d& v) const 
{
    Eigen::Vector4d vh(v.x(), v.y(), v.z(), 1.0);
    return vh.transpose() * Q * vh;
}

Eigen::Vector3d Quadric::findOptimalPosition(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const 
{
    Eigen::Matrix3d A = Q.block<3,3>(0,0);
    Eigen::Vector3d b = -Q.block<3,1>(0,3);
    
    // A is invertible
    Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    if (lu.isInvertible()) 
    {
        return lu.solve(b);
    } 
    else 
    {
        // use midpoint otherwise
        Eigen::Vector3d midpoint = 0.5 * (v1 + v2);
        double error1 = computeError(v1);
        double error2 = computeError(v2);
        double errorMid = computeError(midpoint);
        
        if (error1 <= error2 && error1 <= errorMid) return v1;
        else if (error2 <= errorMid) return v2;
        else return midpoint;
    }
}

bool EdgeCollapse::operator>(const EdgeCollapse& other) const 
{
    return cost > other.cost;
}

Eigen::Vector4d computePlaneEquation(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) 
{
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;
    Eigen::Vector3d normal = v1.cross(v2).normalized();
    double d = -normal.dot(p1);
    return Eigen::Vector4d(normal.x(), normal.y(), normal.z(), d);
}


//
int Simplifier::get_halfedge_between_vertices(Mesh_connectivity& mesh, const int v0, const int v1)
{
	Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh.vertex_ring_at(v0);

	int answer = mesh.invalid_index;

	do
	{

		if(ring_iter.half_edge().origin().index() == v1)
		{
			answer = ring_iter.half_edge().twin().index();
		}
	} while(ring_iter.advance());

	if(answer != mesh.invalid_index)
	{
		assert(mesh.half_edge_at(answer).origin().index() == v0);
		assert(mesh.half_edge_at(answer).dest().index() == v1);
	}

	return answer;
}



bool Simplifier::simplify_once()
{
    return simplify_try_before_commit(mesh().n_active_faces() - 1);
}

bool Simplifier::simplify_test_ahead(int num_faces_to_simplify)
{
    
    int original_num_faces = mesh().n_active_faces();
    // int target_faces = original_num_faces / 2;
    int target_faces  = original_num_faces - num_faces_to_simplify;
    
    // std::cout << "Starting QSlim simplification from " << original_num_faces 
    //           << " faces to " << target_faces << " faces" << std::endl;
    
    std::map<int, Quadric> vertex_quadrics;
    initializeVertexQuadrics(vertex_quadrics);
    
    std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>> edge_queue;
    std::set<std::pair<int,int>> processed_edges; // To avoid duplicate edges
    
    buildInitialEdgeQueue(vertex_quadrics, edge_queue, processed_edges);
    
    // std::cout << "Initial edge queue size: " << edge_queue.size() << std::endl;
    
    // edge collapses
    int collapses_performed = 0;
    while (mesh().n_active_faces() > target_faces && !edge_queue.empty()) 
    {
        EdgeCollapse collapse = edge_queue.top();
        edge_queue.pop();
      
        
        if (!isValidCollapse(collapse)) 
        {
  
          continue;  
        }

        if (performEdgeCollapse(mesh(), collapse, vertex_quadrics)) 
        {
            collapses_performed++;
            
            updateEdgeQueue(collapse.v1, vertex_quadrics, edge_queue, processed_edges);
            
            if (collapses_performed % 1 == 0) 
            {
                std::cout << collapse.v1 << " " << collapse.v2 << std::endl;
                std::cout << "Performed " << collapses_performed << " collapses, "
                          << mesh().n_active_faces() << " faces remaining" << std::endl;
            }
        }
            
    }
    
    std::cout << "QSlim simplification completed. Performed " << collapses_performed 
              << " edge collapses. Final face count: " << mesh().n_active_faces() << std::endl;
    
    return true;
}

// not use anymore
bool Simplifier::simplify_try_before_commit(int target_faces)
{
    initializeProxyMesh();
    
    // Store original number of faces
    int original_num_faces = mesh().n_active_faces();
    // int target_faces = original_num_faces / 2;
    // int target_faces  = 0;
    
    std::cout << "Starting QSlim simplification from " << original_num_faces 
              << " faces to " << target_faces << " faces" << std::endl;
    
    // Step 1: Initialize quadrics for all vertices
    std::map<int, Quadric> vertex_quadrics;
    initializeVertexQuadrics(vertex_quadrics);
    
    // Step 2: Compute initial edge collapse costs and build priority queue
    std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>> edge_queue;
    std::set<std::pair<int,int>> processed_edges; // To avoid duplicate edges
    
    buildInitialEdgeQueue(vertex_quadrics, edge_queue, processed_edges);
    
    std::cout << "Initial edge queue size: " << edge_queue.size() << std::endl;
    
    // Step 3: Perform edge collapses
    int collapses_performed = 0;
    while (mesh().n_active_faces() > target_faces && !edge_queue.empty()) 
    {
        EdgeCollapse collapse = edge_queue.top();
        edge_queue.pop();
        
        // Check if the edge still exists and is valid for collapse
        
        if (!isValidCollapse(collapse)) 
        {
            continue;
        }

        // try the collapse on the proxy mesh
        if (performEdgeCollapse(proxy_mesh(), collapse, vertex_quadrics)) 
        {
            if (checkProxyManifold(proxy_mesh()))
            {
                if (performEdgeCollapse(mesh(), collapse, vertex_quadrics)) 
                {
                    collapses_performed++;
                    
                    updateEdgeQueue(collapse.v1, vertex_quadrics, edge_queue, processed_edges);
                    
                    if (collapses_performed % 1 == 0) 
                    {
                        std::cout << collapse.v1 << " " << collapse.v2 << std::endl;
                        std::cout << "Performed " << collapses_performed << " collapses, "
                                 << mesh().n_active_faces() << " faces remaining" << std::endl;
                    }
                }else{
                    // This should not happen if the proxy collapse succeeded
                    std::cout << "Unexpected: Collapse " << collapse.v1 << " -> " << collapse.v2 
                             << " failed on original mesh after succeeding on proxy" << std::endl;
                    copyMeshToProxy();
                }
            } 
            else 
            {
                std::cout << "Collapse " << collapse.v1 << " -> " << collapse.v2 
                         << " rejected due to manifold violation" << std::endl;
                copyMeshToProxy();
                
            }
        }
        else{
            std::cout << "Collapse " << collapse.v1 << " -> " << collapse.v2 
                      << " could not be performed on proxy mesh" << std::endl;
            copyMeshToProxy();
        }
    }
    
    std::cout << "QSlim simplification completed. Performed " << collapses_performed 
              << " edge collapses. Final face count: " << mesh().n_active_faces() << std::endl;
    
    return true;
}

void Simplifier::initializeVertexQuadrics(std::map<int, Quadric>& vertex_quadrics)
{
    // initialize with zero quadrics
    for (int v_idx = 0; v_idx < mesh().n_total_vertices(); ++v_idx) 
    {
        auto vertex = mesh().vertex_at(v_idx);
        if (vertex.is_active()) 
        {
            vertex_quadrics[v_idx] = Quadric();
        }
    }
    
    // accumulate quadrics from adjacent faces
    for (int f_idx = 0; f_idx < mesh().n_total_faces(); ++f_idx) 
    {
        auto face = mesh().face_at(f_idx);
        if (!face.is_active()) continue;
        
        auto he = face.half_edge();
        int v1_idx = he.origin().index();
        int v2_idx = he.next().origin().index();  
        int v3_idx = he.next().next().origin().index();
        
        Eigen::Vector3d p1 = mesh().vertex_at(v1_idx).xyz();
        Eigen::Vector3d p2 = mesh().vertex_at(v2_idx).xyz();
        Eigen::Vector3d p3 = mesh().vertex_at(v3_idx).xyz();
        
        // get plane equation
        Eigen::Vector4d plane = computePlaneEquation(p1, p2, p3);
        Quadric face_quadric(plane);

        // add to vertex quadrics
        vertex_quadrics[v1_idx] = vertex_quadrics[v1_idx] + face_quadric;
        vertex_quadrics[v2_idx] = vertex_quadrics[v2_idx] + face_quadric;
        vertex_quadrics[v3_idx] = vertex_quadrics[v3_idx] + face_quadric;
    }
}

// Build initial priority queue of edge collapses
void Simplifier::buildInitialEdgeQueue(const std::map<int, Quadric>& vertex_quadrics,
                                     std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>>& edge_queue,
                                     std::set<std::pair<int,int>>& processed_edges)
{
    for (int he_idx = 0; he_idx < mesh().n_total_half_edges(); ++he_idx) 
    {
        auto he = mesh().half_edge_at(he_idx);
        if (!he.is_active()) continue;
        
        int v1 = he.origin().index();
        int v2 = he.dest().index();
        
        if (v1 > v2) std::swap(v1, v2);
        
        if (processed_edges.find({v1, v2}) != processed_edges.end()) 
        {
            continue;
        }
        processed_edges.insert({v1, v2});
        
        EdgeCollapse collapse = computeEdgeCollapse(v1, v2, he_idx, vertex_quadrics);
        if (collapse.cost < std::numeric_limits<double>::infinity()) 
        {
            edge_queue.push(collapse);
        }
    }
}

//  cost and optimal position for collapsing an edge
EdgeCollapse Simplifier::computeEdgeCollapse(int v1, int v2, int he_idx, const std::map<int, Quadric>& vertex_quadrics)
{
    EdgeCollapse collapse;
    collapse.he_index = he_idx;
    collapse.v1 = v1;
    collapse.v2 = v2;
    
    auto it1 = vertex_quadrics.find(v1);
    auto it2 = vertex_quadrics.find(v2);
    
    if (it1 == vertex_quadrics.end() || it2 == vertex_quadrics.end()) 
    {
        collapse.cost = std::numeric_limits<double>::infinity();
        return collapse;
    }
    
    Quadric combined = it1->second + it2->second;
    
    // optimal position
    Eigen::Vector3d pos1 = mesh().vertex_at(v1).xyz();
    Eigen::Vector3d pos2 = mesh().vertex_at(v2).xyz();
    
    collapse.new_pos = combined.findOptimalPosition(pos1, pos2);
    collapse.cost = combined.computeError(collapse.new_pos);
    
    return collapse;
}

// if an edge collapse is still valid
bool Simplifier::isValidCollapse(const EdgeCollapse& collapse)
{
    auto v1 = mesh().vertex_at(collapse.v1);
    auto v2 = mesh().vertex_at(collapse.v2);
    
    if (!v1.is_active() || !v2.is_active()) 
    {
        return false;
    }
    
    // whether halfedge still exists
    if (collapse.he_index >= 0 && collapse.he_index < mesh().n_total_half_edges()) 
    {
        auto he = mesh().half_edge_at(collapse.he_index);
        if (he.is_active() && 
            ((he.origin().index() == collapse.v1 && he.dest().index() == collapse.v2) ||
             (he.origin().index() == collapse.v2 && he.dest().index() == collapse.v1))) 
        {
            // two-neighbor for the edge
            return !wouldCreateNonManifold(collapse.v1, collapse.v2);
        }
    }
    
    int he_between = get_halfedge_between_vertices(mesh(), collapse.v1, collapse.v2);
    if (he_between != mesh().invalid_index) 
    {
        return !wouldCreateNonManifold(collapse.v1, collapse.v2);
    }
    
    return false;
}

// perform edge collapse
bool Simplifier::performEdgeCollapse(Mesh_connectivity& mesh, EdgeCollapse& collapse, std::map<int, Quadric>& vertex_quadrics)
{

    int he_idx = get_halfedge_between_vertices(mesh, collapse.v1, collapse.v2);
    if (he_idx == mesh.invalid_index) 
    {
        he_idx = get_halfedge_between_vertices(mesh, collapse.v2, collapse.v1);
    }
    
    if (he_idx == mesh.invalid_index) 
    {
        return false;
    }
    
    auto he = mesh.half_edge_at(he_idx);
    
    // sort to avoid collision
    if (he.origin().index() != collapse.v1) 
    {
        std::swap(collapse.v1, collapse.v2);
    }
    
    // update position
    auto v1 = mesh.vertex_at(collapse.v1);
    v1.data().xyz = collapse.new_pos;
    
    vertex_quadrics[collapse.v1] = vertex_quadrics[collapse.v1] + vertex_quadrics[collapse.v2];
    vertex_quadrics.erase(collapse.v2);
    
    // important
    return collapseEdgeTopology(mesh, he_idx, collapse.v1, collapse.v2);
}


void Simplifier::checkVertexAssociatedHalfEdges(Mesh_connectivity& mesh,int halfedge_idx)
{

    auto halfedge = mesh.half_edge_at(halfedge_idx);
    for (int i = 0; i < mesh.n_total_vertices(); ++i) 
    {
        auto vertex = mesh.vertex_at(i);
        if (vertex.is_active() && vertex.half_edge().index() == halfedge_idx)
        {
            do{
                mesh.vertex_at(i).data().half_edge = vertex.half_edge().prev().twin().index();
                }while(!mesh.half_edge_at(vertex.data().half_edge).is_active());
        }
    }
}


// update connectivity
bool Simplifier::collapseEdgeTopology(Mesh_connectivity & mesh, int he_idx, int v_keep, int v_remove)
{
    auto he = mesh.half_edge_at(he_idx);
    auto he_twin = he.twin();
    
    auto f1 = he.face();
    auto f2 = he_twin.face();
    
    auto he1_next = he.next();
    auto he1_prev = he.prev();
    auto he2_next = he_twin.next();
    auto he2_prev = he_twin.prev();
    
    // /\
    // \/
    auto he1_next_twin = he1_next.twin();
    auto he1_prev_twin = he1_prev.twin();
    auto he2_next_twin = he2_next.twin();
    auto he2_prev_twin = he2_prev.twin();
    
    // halfedges_to_redirect
    auto v_remove_iter = mesh.vertex_at(v_remove);
    auto ring_iter = mesh.vertex_ring_at(v_remove);
    
    std::vector<int> halfedges_to_redirect;
    do 
    {
        int current_he = ring_iter.half_edge().index();
        if (current_he != he.index() && current_he != he_twin.index() &&
            current_he != he1_next.index() && current_he != he1_prev.index() &&
            current_he != he2_next.index() && current_he != he2_prev.index() &&
            current_he != he1_next_twin.index() && current_he != he1_prev_twin.index() &&
            current_he != he2_next_twin.index() && current_he != he2_prev_twin.index())
        {
            halfedges_to_redirect.push_back(current_he);
        }
    } while (ring_iter.advance());


    for (int redirect_he : halfedges_to_redirect) 
    {
        auto redirect_he_iter = mesh.half_edge_at(redirect_he);
        
        redirect_he_iter.twin().data().origin = v_keep;
      }


    // for (int redirect_he : halfedges_to_redirect) 
    // {
    //     auto redirect_he_iter = mesh().half_edge_at(redirect_he);
    //     redirect_he_iter.data().origin = v_keep;
    // }
    
    auto v_keep_ring = mesh.vertex_ring_at(v_keep);
    bool found_valid_he = false;

    if (he1_prev_twin.is_active() && he1_next_twin.is_active()) {
        he1_prev_twin.data().twin = he1_next_twin.index();
        he1_next_twin.data().twin = he1_prev_twin.index();
    }
    
    // For triangle f2: connect the twins of he2_prev and he2_next  
    if (he2_prev_twin.is_active() && he2_next_twin.is_active()) {
        he2_prev_twin.data().twin = he2_next_twin.index();
        he2_next_twin.data().twin = he2_prev_twin.index();
    }


    he2_prev.twin().data().origin = v_keep;

    he.deactivate();
    // find a valid halfedge to associate with v_keep
    checkVertexAssociatedHalfEdges(mesh, he.index());

    he_twin.deactivate();
    checkVertexAssociatedHalfEdges(mesh, he_twin.index());

    he1_next.deactivate();
    checkVertexAssociatedHalfEdges(mesh, he1_next.index());

    he1_prev.deactivate();
    checkVertexAssociatedHalfEdges(mesh, he1_prev.index());

    he2_next.deactivate();
    checkVertexAssociatedHalfEdges(mesh, he2_next.index());

    he2_prev.deactivate();
    checkVertexAssociatedHalfEdges(mesh, he2_prev.index());

    // Deactivate the faces
    if (f1.is_active() && !f1.is_equal(mesh.hole())) 
    {
        f1.deactivate();
    }
    if (f2.is_active() && !f2.is_equal(mesh.hole())) 
    {
        f2.deactivate();
    }

    
    v_remove_iter.deactivate();
    
    return true;
}

void Simplifier::updateEdgeQueue(int vertex_idx, const std::map<int, Quadric>& vertex_quadrics,
                                std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>>& edge_queue,
                                std::set<std::pair<int,int>>& processed_edges)
{
    auto ring_iter = mesh().vertex_ring_at(vertex_idx);
    
    std::set<int> neighbors;
    do 
    {
        int neighbor = ring_iter.half_edge().origin().index();
        if (neighbors.find(neighbor) == neighbors.end()) 
        {
                       neighbors.insert(neighbor);
            
            int v1 = vertex_idx;
            int v2 = neighbor;
            if (v1 > v2) std::swap(v1, v2);
            
            if (processed_edges.find({v1, v2}) == processed_edges.end()) 
            {
                processed_edges.insert({v1, v2});
                EdgeCollapse collapse = computeEdgeCollapse(v1, v2, ring_iter.half_edge().index(), vertex_quadrics);
                if (collapse.cost < std::numeric_limits<double>::infinity()) 
                {
                    edge_queue.push(collapse);
                }
            }
        }
    } while (ring_iter.advance());
}

bool Simplifier::wouldCreateNonManifold(int v1, int v2)
{
    auto vertex1 = mesh().vertex_at(v1);
    auto vertex2 = mesh().vertex_at(v2);
    
    if (!vertex1.is_active() || !vertex2.is_active()) {
        std::cout<<"One of the vertices is inactive: " << v1 << " or " << v2 << std::endl;
      return true; // Can't collapse inactive vertices
    }
    

    // important: only allow 2 vectex neighbors for the edge
    std::set<int> neighbors_v1, neighbors_v2, shared_neighbors;
    
    auto ring1 = mesh().vertex_ring_at(v1);
    do {
        int neighbor = ring1.half_edge().origin().index();
        if (neighbor != v2) { // Don't include v2 as a neighbor of v1
            neighbors_v1.insert(neighbor);
        }
    } while (ring1.advance());
    
    auto ring2 = mesh().vertex_ring_at(v2);
    do {
        int neighbor = ring2.half_edge().origin().index();
        if (neighbor != v1) { // Don't include v1 as a neighbor of v2
            neighbors_v2.insert(neighbor);
        }
    } while (ring2.advance());
    
    std::set_intersection(neighbors_v1.begin(), neighbors_v1.end(),
                         neighbors_v2.begin(), neighbors_v2.end(),
                         std::inserter(shared_neighbors, shared_neighbors.begin()));
    
    // Basic manifold checks:
    
    // if (neighbors_v1.size() < 2 || neighbors_v2.size() < 2) {
    //     return true; // Vertices with degree < 2 are problematic
    // }
    
    if (shared_neighbors.size() != 2) {
      std::cout << "Non-manifold detected: shared neighbors > 2 between vertices " << v1 << " and " << v2 << std::endl;
        
      return true;
    }
    
   
    return false; 
}

// not use anymore, initialize proxy mesh as a copy
void Simplifier::initializeProxyMesh()
{
    try {
        _proxy_mesh.copy(mesh());
        std::cout << "Proxy mesh initialized successfully with " 
                  << _proxy_mesh.n_active_vertices() << " vertices and "
                  << _proxy_mesh.n_active_faces() << " faces" << std::endl;
    }
    catch (...) {
        std::cerr << "Failed to copy mesh to proxy mesh, trying alternative approach" << std::endl;
        
        // Alternative approach: use build_from_triangles
        try {
            std::vector<double> vertices;
            std::vector<int> triangles;
            std::map<int, int> vertex_map; // old_index -> new_index
            
            int new_vertex_count = 0;
            for (int v_idx = 0; v_idx < mesh().n_total_vertices(); ++v_idx) 
            {
                auto vertex = mesh().vertex_at(v_idx);
                if (vertex.is_active()) 
                {
                    vertex_map[v_idx] = new_vertex_count++;
                    auto pos = vertex.xyz();
                    vertices.push_back(pos.x());
                    vertices.push_back(pos.y());
                    vertices.push_back(pos.z());
                }
            }
            
            for (int f_idx = 0; f_idx < mesh().n_total_faces(); ++f_idx) 
            {
                auto face = mesh().face_at(f_idx);
                if (face.is_active()) 
                {
                    auto he = face.half_edge();
                    int v1_orig = he.origin().index();
                    int v2_orig = he.next().origin().index();
                    int v3_orig = he.next().next().origin().index();
                    
                    if (vertex_map.find(v1_orig) != vertex_map.end() &&
                        vertex_map.find(v2_orig) != vertex_map.end() &&
                        vertex_map.find(v3_orig) != vertex_map.end()) 
                    {
                        triangles.push_back(vertex_map[v1_orig]);
                        triangles.push_back(vertex_map[v2_orig]);
                        triangles.push_back(vertex_map[v3_orig]);
                    }
                }
            }
            
            // Build proxy mesh from extracted data
            _proxy_mesh.clear();
            _proxy_mesh.build_from_triangles(vertices, triangles);
            
            std::cout << "Proxy mesh built from triangles with " 
                      << _proxy_mesh.n_active_vertices() << " vertices and "
                      << _proxy_mesh.n_active_faces() << " faces" << std::endl;
        }
        catch (...) {
            std::cerr << "Failed to build proxy mesh from triangles" << std::endl;
            return;
        }
    }
}

void Simplifier::copyMeshToProxy()
{
    initializeProxyMesh();
}

void Simplifier::copyProxyToMesh()
{
    
    for (int v_idx = 0; v_idx < std::min(mesh().n_total_vertices(), _proxy_mesh.n_total_vertices()); ++v_idx) 
    {
        auto orig_vertex = mesh().vertex_at(v_idx);
        auto proxy_vertex = _proxy_mesh.vertex_at(v_idx);
        
        if (orig_vertex.is_active() && proxy_vertex.is_active()) 
        {
            orig_vertex.data().xyz = proxy_vertex.xyz();
        }
    }
}

// not use anymore
bool Simplifier::checkProxyManifold(Mesh_connectivity& mesh)
{
    
    // Check edge-face relationships
    std::map<std::pair<int,int>, int> edge_face_count;
    
    for (int f_idx = 0; f_idx < mesh.n_total_faces(); ++f_idx) 
    {
        auto face = mesh.face_at(f_idx);
        if (!face.is_active()) continue;
        
        // Get the three vertices of the triangle
        auto he = face.half_edge();
        int v1 = he.origin().index();
        int v2 = he.next().origin().index();
        int v3 = he.next().next().origin().index();
        
        // Check each edge of the triangle
        std::vector<std::pair<int,int>> edges = {
            {std::min(v1,v2), std::max(v1,v2)},
            {std::min(v2,v3), std::max(v2,v3)},
            {std::min(v3,v1), std::max(v3,v1)}
        };
        
        for (auto& edge : edges) 
        {
            edge_face_count[edge]++;
            
            // If any edge is shared by more than 2 faces, it's non-manifold
            if (edge_face_count[edge] > 2) 
            {
                std::cout << "Non-manifold detected: Edge (" << edge.first << "," << edge.second 
                         << ") shared by " << edge_face_count[edge] << " faces" << std::endl;
                return false;
            }
        }
    }

//     for(int i=0; i < mesh.n_total_half_edges(); ++i)
//     {
//         auto he = mesh.half_edge_at(i);
//         if(!he.is_active())
//             continue;

//         auto he_twin = he.twin();
//         if(!he_twin.is_active())
//             continue;

//         edge_idx = {std::min(he.origin().index(), he.dest().index()), std::max(he.origin().index(), he.dest().index())};
//         edge_face_count[edge_idx]++;
        
//         if(edge_face_count[edge_idx] > 2)
//         {
//             std::cout << "Non-manifold detected: Half-edge (" << he.index() << "," << he_twin.index() 
//                       << ") shared by more than 2 faces" << std::endl;
//             return false;
//         }

//     }
// }

    
    return true; // All checks passed - mesh appears manifold
}

// Get the edges with the lowest collapse cost
std::vector<EdgeCollapse> Simplifier::getLowestCostEdges(int num_edges)
{
    std::vector<EdgeCollapse> result;
    
    // initialize quadrics for all vertices
    std::map<int, Quadric> vertex_quadrics;
    initializeVertexQuadrics(vertex_quadrics);
    
    // priority queue 
    std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>> edge_queue;
    std::set<std::pair<int,int>> processed_edges;
    
    buildInitialEdgeQueue(vertex_quadrics, edge_queue, processed_edges);
    
    // Extract the lowest cost edges
    for (int i = 0; i < num_edges && !edge_queue.empty(); ++i) {
        EdgeCollapse collapse = edge_queue.top();
        edge_queue.pop();
        
  
        if (isValidCollapse(collapse)) {
            result.push_back(collapse);
        } else {
            --i;
        }
    }
    
    return result;
}

Eigen::Matrix4Xf Simplifier::createEdgeCostVertexColors(const std::vector<EdgeCollapse>& lowest_cost_edges,
                                                        Mesh_connectivity::Defragmentation_maps& defrag)
{
    Eigen::Matrix4Xf colors(4, mesh().n_active_vertices());
    
    // Default color: light gray
    for (int i = 0; i < mesh().n_active_vertices(); ++i) {
        colors(0, i) = 0.8f; 
        colors(1, i) = 0.8f; 
        colors(2, i) = 0.8f; 
        colors(3, i) = 1.0f; 
    }
    
    // Priority tracking: 0 = default, 1 = yellow, 2 = orange, 3 = red
    std::vector<int> vertex_priority(mesh().n_active_vertices(), 0);
    
    // Color vertices of lowest cost edges
    std::vector<Eigen::Vector3f> edge_colors = {
        {1.0f, 0.0f, 0.0f}, // Red for lowest cost
        {1.0f, 0.5f, 0.0f}, // Orange for second lowest
        {1.0f, 1.0f, 0.0f}  // Yellow for third lowest
    };
    
    // Process in reverse order (yellow first) so red can overwrite
    for (int i = lowest_cost_edges.size() - 1; i >= 0 && i < (int)edge_colors.size(); --i) {
        const EdgeCollapse& collapse = lowest_cost_edges[i];
        int priority = 3 - i; // Red=3, Orange=2, Yellow=1
        
        // Color both vertices of the edge if priority is higher
        int v1_continuous = defrag.old2new_vertices[collapse.v1];
        int v2_continuous = defrag.old2new_vertices[collapse.v2];
        
        if (v1_continuous >= 0 && v1_continuous < mesh().n_active_vertices()) {
            if (vertex_priority[v1_continuous] < priority) {
                colors(0, v1_continuous) = edge_colors[i](0);
                colors(1, v1_continuous) = edge_colors[i](1);
                colors(2, v1_continuous) = edge_colors[i](2);
                vertex_priority[v1_continuous] = priority;
            }
        }
        
        if (v2_continuous >= 0 && v2_continuous < mesh().n_active_vertices()) {
            if (vertex_priority[v2_continuous] < priority) {
                colors(0, v2_continuous) = edge_colors[i](0);
                colors(1, v2_continuous) = edge_colors[i](1);
                colors(2, v2_continuous) = edge_colors[i](2);
                vertex_priority[v2_continuous] = priority;
            }
        }
    }
    
    return colors;
}

// Create face colors to highlight faces adjacent to lowest cost edges
Eigen::Matrix4Xf Simplifier::createEdgeCostFaceColors(const std::vector<EdgeCollapse>& lowest_cost_edges,
                                                      Mesh_connectivity::Defragmentation_maps& defrag)
{
    Eigen::Matrix4Xf colors(4, mesh().n_active_faces());
    
    // Default color: light gray
    for (int i = 0; i < mesh().n_active_faces(); ++i) {
        colors(0, i) = 0.9f; // R
        colors(1, i) = 0.9f; // G
        colors(2, i) = 0.9f; // B
        colors(3, i) = 1.0f; // A
    }
    
    // Priority tracking: 0 = default, 1 = yellow, 2 = orange, 3 = red
    std::vector<int> face_priority(mesh().n_active_faces(), 0);
    
    // Color faces adjacent to lowest cost edges
    std::vector<Eigen::Vector3f> edge_colors = {
        {1.0f, 0.0f, 0.0f}, // Red for lowest cost
        {1.0f, 0.5f, 0.0f}, // Orange for second lowest  
        {1.0f, 1.0f, 0.0f}  // Yellow for third lowest
    };
    
    // Process in reverse order (yellow first) so red can overwrite
    for (int i = lowest_cost_edges.size() - 1; i >= 0 && i < (int)edge_colors.size(); --i) {
        const EdgeCollapse& collapse = lowest_cost_edges[i];
        int priority = 3 - i; // Red=3, Orange=2, Yellow=1
        
        // Find faces that contain this edge
        if (collapse.he_index >= 0 && collapse.he_index < mesh().n_total_half_edges()) {
            auto he = mesh().half_edge_at(collapse.he_index);
            if (he.is_active()) {
                // Color the face containing this half-edge
                if (!he.face().is_equal(mesh().hole())) {
                    int face_idx = he.face().index();
                    int face_continuous = defrag.old2new_faces[face_idx];
                    
                    if (face_continuous >= 0 && face_continuous < mesh().n_active_faces()) {
                        if (face_priority[face_continuous] < priority) {
                            colors(0, face_continuous) = edge_colors[i](0);
                            colors(1, face_continuous) = edge_colors[i](1);
                            colors(2, face_continuous) = edge_colors[i](2);
                            face_priority[face_continuous] = priority;
                        }
                    }
                }
                
                // Color the face containing the twin half-edge
                auto he_twin = he.twin();
                if (he_twin.is_active() && !he_twin.face().is_equal(mesh().hole())) {
                    int face_idx = he_twin.face().index();
                    int face_continuous = defrag.old2new_faces[face_idx];
                    
                    if (face_continuous >= 0 && face_continuous < mesh().n_active_faces()) {
                        if (face_priority[face_continuous] < priority) {
                            colors(0, face_continuous) = edge_colors[i](0);
                            colors(1, face_continuous) = edge_colors[i](1);
                            colors(2, face_continuous) = edge_colors[i](2);
                            face_priority[face_continuous] = priority;
                        }
                    }
                }
            }
        }
    }
    
    return colors;
}

} // end of mohecore
} // end of minimesh