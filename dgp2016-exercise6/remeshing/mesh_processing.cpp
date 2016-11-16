//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#include "mesh_processing.h"
#include <set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using std::min;
using std::max;
using std::cout;
using std::endl;
using namespace Eigen;

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

MeshProcessing::~MeshProcessing() {
    // TODO
}

void MeshProcessing::remesh (const REMESHING_TYPE &remeshing_type,
                             const int &num_iterations) {
    calc_weights ();
    calc_mean_curvature ();
    calc_uniform_mean_curvature ();
    calc_gauss_curvature ();
    calc_target_length (remeshing_type);

    // main remeshing loop
    for (int i = 0; i < num_iterations; ++i){
        split_long_edges ();
        collapse_short_edges ();
        equalize_valences ();
        tangential_relaxation ();
        std::cout << "remesh number " << i << " done" << std::endl;
    }
}

void MeshProcessing::calc_target_length (const REMESHING_TYPE &remeshing_type){
    Scalar                   length;
    Scalar                   mean_length;
    Scalar                   user_specified_target_length;

    Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0);
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    Mesh::Vertex_property<Scalar> target_new_length  = mesh_.vertex_property<Scalar>("v:newlength", 0);

    // compute mean edge length for future use
    mean_length = 0.0;
    int n = 0;
    for(auto e: mesh_.edges()){
        if(e.is_valid()){
            mean_length += mesh_.edge_length(e);
            n++;
        }
    }
    mean_length /= n;

    // Scaled average edge length
    if (remeshing_type == AVERAGE){
        // assign target length for each vertex v
        for(auto v:mesh_.vertices()){
            target_length[v] = mean_length;
        }
    }

    // Curvature-based adaptive remeshing
    else if (remeshing_type == CURV){
        // instantiate min curvature to 0.0
        double min_curv = 0.0;
        for(auto v:mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                // find min curvature
                if(curvature[v] < min_curv){
                    min_curv = curvature[v];
                }
            }
        }

        // assign target length for each vertex v
        for(auto v:mesh_.vertices()){
            length = 1.0;
            if(!mesh_.is_boundary(v)){
                // get curvature of current vertex
                length = curvature[v];
            }
            // update target length with the curvature of current vertex + min value (to avoid negative target length)
            target_length[v] = length+(std::abs(min_curv));
        }


        // smooth desired length (using uniform smooth)
        for (int i = 0; i < 5; i++) {
            for(auto v: mesh_.vertices()){
                if(!mesh_.is_boundary(v)){
                    unsigned int N = mesh_.valence(v);
                    double sum = 0.0;
                    // iterate over all neighbors of v
                    for(auto neighbor: mesh_.vertices(v)){
                        sum += target_length[neighbor];
                    }
                    double update = (sum/(double)N) - target_length[v];
                    // compute the new target length
                    target_length[v] = target_length[v] + 0.5 * update;
                }
            }
        }

        // rescale desired length:
        double sum = 0.0;
        int n = 0;
        for(auto v: mesh_.vertices()){
            sum += target_length[v];
            n++;
        }

        user_specified_target_length = mean_length;

        // rescale the target length of each vertex so that the mean of the new target lengths equals the user specified target length
        for(auto v: mesh_.vertices()){
            // compute new target length
            target_new_length[v] = n*user_specified_target_length*target_length[v]/sum;
            target_length[v] = target_new_length[v];
        }
    }

    // Height-based
    else if (remeshing_type == HEIGHT){
        // instantiate min height to 0.0
        double min_height = 0.0;
        for(auto v:mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                // find min height
                if(mesh_.position(v)[1] < min_height){
                    min_height = mesh_.position(v)[1];
                }
            }
        }

        // assign target length for each vertex v
        for(auto v:mesh_.vertices()){
            // get height of current vertex
            length = mesh_.position(v)[1];

            // update target length with the height of current vertex + min value (to avoid negative target length)
            target_length[v] = length+(std::abs(min_height));
        }

        // rescale desired length:
        double sum = 0.0;
        int n = 0;
        for(auto v: mesh_.vertices()){
            sum += target_length[v];
            n++;
        }

        user_specified_target_length = 10.0;

        // rescale the target length of each vertex so that the mean of the new target lengths equals the user specified target length
        for(auto v: mesh_.vertices()){
            // compute new target length
            target_new_length[v] = n*user_specified_target_length*target_length[v]/sum;
            target_length[v] = target_new_length[v] + 1.0;
        }

    }
    std::cout << "Calc target length done" << std::endl;
}

void MeshProcessing::split_long_edges (){
    Mesh::Vertex   v0, v1, v;
    bool            finished;
    int             i;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished=false, i=0; !finished && i<100; ++i){
        // set to true at the beginning and later set to false if at least one edge is split during this iteration
        finished = true;
        // iterate over each edge e
        for (auto e:mesh_.edges()){
            // find the two endpoints of e (v1, v2)
            v0 = mesh_.vertex(e, 0);
            v1 = mesh_.vertex(e, 1);
            // compute e's edge target length
            double target_length_e = ((target_length[v0] + target_length[v1])/2.0);
            //split e if its length is bigger than 4/3 times its target length
            if(mesh_.edge_length(e) > target_length_e*(4.0/3.0)){
                finished = false;
                //split e and add a new vertex v between v0 and v1
                v = mesh_.add_vertex(Point((mesh_.position(v0)+mesh_.position(v1))/2.0));
                mesh_.split(e, v);

                //compute and store v's normal
                normals[v] = mesh_.compute_vertex_normal(v);
                // interpolate target edge length of new vertex v
                target_length[v] = (target_length[v0]+target_length[v1])/2.0;
            }
        }
    }
    std::cout << "Split long edges done" << std::endl;
}
void MeshProcessing::collapse_short_edges (){
    Mesh::Vertex   v0, v1;
    Mesh::Halfedge  h01, h10;
    bool            finished, b0, b1;
    int             i;
    bool            hcol01, hcol10;

    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    for (finished=false, i=0; !finished && i<100; ++i){
        // finished = true if we don't collapse any edge
        // finished = false in the other case
        finished = true;

        for (auto e: mesh_.edges()){
            // might already be deleted
            if (!mesh_.is_deleted(e)){
                //we find the two end vertices of the edge
                v0 = mesh_.vertex(e,0);
                v1 = mesh_.vertex(e,1);
                // compute e's target length from v0's and v1's.
                double target_length_e = ((target_length[v0] + target_length[v1])/2.0);

                // check if the edge is smaller the (4/5) of the edge target length
                if (mesh_.edge_length(e) < (4.0/5.0)*target_length_e){
                    // we compute every useful variable corresponding to the edge
                    b0 = mesh_.is_boundary(v0);
                    b1 = mesh_.is_boundary(v1);
                    h01 = mesh_.find_halfedge(v0,v1);
                    h10 = mesh_.find_halfedge(v1,v0);
                    hcol01 = mesh_.is_collapse_ok(h01);
                    hcol10 = mesh_.is_collapse_ok(h10);

                    // If neither endpoint is a boundary vertex or if both endpoints are boundary vertices, both halfedges are potentially collapsable
                    if ((!b0 && !b1) || (b0 && b1)) {
                         //If both halfedges are collapsable, collapse the lower valence vertex into the higher one
                        if (hcol01 && hcol10){
                            finished = false;
                            if (mesh_.valence(v0) < mesh_.valence(v1)){
                                mesh_.collapse(h01);
                            } else {
                                mesh_.collapse(h10);
                            }
                        // If exactly one of the halfedges is collapsable, collapse only this halfedge
                        } else if (hcol01) {
                            mesh_.collapse(h01);
                            finished = false;
                        } else if (hcol10){
                            mesh_.collapse(h10);
                            finished = false;
                        }
                    // if exactly one endpoint of e is a boundary v, we only collapse the halfedge going from the non-b. to the b. vertex
                    } else if (!b0){
                        if (hcol01){
                            mesh_.collapse(h01);
                            finished = false;
                        }
                    } else if (!b1){
                        if (hcol10){
                            mesh_.collapse(h10);
                            finished = false;
                        }
                    }
                }
            }
        }
    }
    // clean the deleted edges/vertices/faces
    mesh_.garbage_collection();

    if (i==100) std::cerr << "collapse break\n";
    std::cout << "Collapse short edges done" << std::endl;
}

void MeshProcessing::equalize_valences (){
    Mesh::Vertex   v0, v1, v2, v3;
    Mesh::Halfedge   h;
    int             val0, val1, val2, val3;
    int             val_opt0, val_opt1, val_opt2, val_opt3;
    int             ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool            finished;
    int             i;


    // flip all edges
    for (finished=false, i=0; !finished && i<100; ++i){

        // finished = true if we don't flip any edge
        // finished = false in the other case
        finished = true;

        // iterate on edges of the mesh
        for (auto e: mesh_.edges()){
            // check if current edge is boundary
            if (!mesh_.is_boundary(e))
            {
                // find the two end vertices of the edge
                v0 = mesh_.vertex(e,0);
                v1 = mesh_.vertex(e,1);

                // get the halfedge between v0 and v1
                h = mesh_.find_halfedge(v0,v1);

                // get v2 and v3 using the next halfedges from the h and its opposite
                v2 = mesh_.to_vertex(mesh_.next_halfedge(h));
                v3 = mesh_.to_vertex(mesh_.next_halfedge(mesh_.opposite_halfedge(h)));

                // get valences of each vertex
                val0 = mesh_.valence(v0);
                val1 = mesh_.valence(v1);
                val2 = mesh_.valence(v2);
                val3 = mesh_.valence(v3);

                // check what is the optimal valence for each vertex
                (mesh_.is_boundary(v0)) ? val_opt0 = 4 : val_opt0 = 6;
                (mesh_.is_boundary(v1)) ? val_opt1 = 4 : val_opt1 = 6;
                (mesh_.is_boundary(v2)) ? val_opt2 = 4 : val_opt2 = 6;
                (mesh_.is_boundary(v3)) ? val_opt3 = 4 : val_opt3 = 6;

                // compute the valence difference
                ve0 = val0 - val_opt0;
                ve1 = val1 - val_opt1;
                ve2 = val2 - val_opt2;
                ve3 = val3 - val_opt3;

                // compute the squarred sum of these differences
                ve_before = pow(ve0,2) + pow(ve1,2) + pow(ve2,2) + pow(ve3,2);

                // get valences of each vertex after simulating the flip
                val0 -= 1;
                val1 -= 1;
                val2 += 1;
                val3 += 1;

                // compute the valence difference
                ve0 = val0 - val_opt0;
                ve1 = val1 - val_opt1;
                ve2 = val2 - val_opt2;
                ve3 = val3 - val_opt3;

                // compute the squarred sum of these differences
                ve_after = pow(ve0,2) + pow(ve1,2) + pow(ve2,2) + pow(ve3,2);

                // compare sums and flip if better
                if(ve_after < ve_before){
                    if(mesh_.is_flip_ok(e)){
                        // flip current edge
                        mesh_.flip(e);
                        finished = false;
                    }
                }
            }
        }
    }
    if (i==100) std::cerr << "flip break\n";
    std::cout << "Equalize valence done" << std::endl;
}

void MeshProcessing::tangential_relaxation (){
    int    valence;
    Point     u, n;
    Point     laplace;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");

    // smooth
    for (int iters=0; iters<10; ++iters){

        // iterate on all vertices in the mesh
        for (auto v: mesh_.vertices()){
            // check if current vertex is boundary
            if (!mesh_.is_boundary(v)){

                // get the normal of current vertex
                n = normals[v];

                // Find its valence
                valence = mesh_.valence(v);

                Point sum(0.0, 0.0, 0.0);
                // iterate through neighbors to get the sum of their positions
                for(auto vn:mesh_.vertices(v)){
                    sum += mesh_.position(vn);
                }

                // compute the laplace point
                laplace = (1/(double) valence) * sum;

                // instantiate the update vector
                Point update_vector(0.0, 0.0, 0.0);

                // compute the update vector
                update_vector[0] = (((1-pow(n[0],2)) * (laplace[0] - mesh_.position(v)[0])) - ((n[0]*n[1])*(laplace[1] - mesh_.position(v)[1])) - ((n[0]*n[2])*(laplace[2] - mesh_.position(v)[2])));
                update_vector[1] = ((-(n[0]*n[1]) * (laplace[0] - mesh_.position(v)[0])) + ((1-pow(n[1],2))*(laplace[1] - mesh_.position(v)[1])) - ((n[1]*n[2])*(laplace[2] - mesh_.position(v)[2])));
                update_vector[2] = ((-(n[0]*n[2]) * (laplace[0] - mesh_.position(v)[0])) - ((n[1]*n[2])*(laplace[1] - mesh_.position(v)[1])) + ((1-pow(n[2],2))*(laplace[2] - mesh_.position(v)[2])));

                // get a lambda
                double lambda = 1;

                // compute the update vector
                u = lambda * update_vector;

                // add it in the update property
                update[v] = u;
            }
        }

        // iterate on all vertices in the mesh
        for (auto v: mesh_.vertices()){
            if (!mesh_.is_boundary(v)){
                // refresh its position using the update property
                mesh_.position(v) += update[v];
            }
        }
    }

    std::cout << "Tangential relaxation done" << std::endl;
}

// ========================================================================
// EXERCISE 1.1
// ========================================================================
void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate mean curvature using
    // the length of the uniform Laplacian approximation
    // Save your approximation in unicurvature vertex property of the mesh.
    // ------------- IMPLEMENT HERE ---------
    for(auto v: mesh_.vertices()){
        if(!mesh_.is_boundary(v)){
            unsigned int N = mesh_.valence(v);
            Point sum(0.0, 0.0, 0.0);
            // iterate over all neighbors of v
            for(auto neighbor: mesh_.vertices(v)){
                sum += (mesh_.position(neighbor) - mesh_.position(v));
            }
            v_unicurvature[v] = (norm(sum)/(double)N);
        }
    }
}

// ========================================================================
// EXERCISE 1.2
// ========================================================================
void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // ------------- IMPLEMENT HERE ---------
    // run calc_weights function to fill e_weight and v_weight arrays
    calc_weights();
    //iterate over each vertex v
    for(auto v: mesh_.vertices()){
        if(!mesh_.is_boundary(v)){
            Point sum(0.0, 0.0, 0.0);
            double w = v_weight[v];
            // iterate over each neighbor vi
            for(auto vi: mesh_.vertices(v)){
                // get the edge between v and vi and get its weight from e_weight
                auto ei = mesh_.find_edge(v, vi);
                double wi = e_weight[ei];
                sum += wi*(mesh_.position(vi)-mesh_.position(v));
            }
            sum *= w;
            v_curvature[v] = norm(sum);
        }
    }
}

// ========================================================================
// EXERCISE 1.3
// ========================================================================
void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // ------------- IMPLEMENT HERE ---------
    // run calc_verices_weights function to fill v_weight arrays
    calc_vertices_weights();
    // Iterate over each vertex v

    for(auto v: mesh_.vertices()){
        if(!mesh_.is_boundary(v)){
            double angle_sum = 0;
            double w = v_weight[v];
            // as w = 1/2A then A = 1/2w
            double area = 2*w;
            // iterate over all neighbors of v
            for(auto v1: mesh_.vertices(v)){
                // iterate over all neighbors of v1
                for(auto v2: mesh_.vertices(v1)){
                    // if there exists an edge between a neighbor of v1 and v then we have to take this angle
                    if(mesh_.find_edge(v,v2).is_valid()==1){
                        // create vectors from points
                        auto vec1 = mesh_.position(v1)-mesh_.position(v);
                        auto vec2 = mesh_.position(v2)-mesh_.position(v);
                        // compute the dot product between two adjacent vectors
                        double result = dot(vec1,vec2);
                        // compute acos(dot(A,B) / (norm(A) * norm(B)))
                        angle_sum += acos(double(result/double(norm(vec1)*norm(vec2))));                    }
                }
            }
            // divide the sum of angle by 2 because each angle is added twice in our sum
            angle_sum = angle_sum/2;
            // compute curvature using formula given
            v_gauss_curvature[v] = (2.0 * M_PI - angle_sum)/area;
        }
    }
}


void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges())
    {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0))
        {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::compute_mesh_properties() {
    Mesh::Vertex_property<Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
                 8 /* max */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for(auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
                mesh_.position(v).y,
                mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
                vertex_normal[v].y,
                vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
                v_color_valence[v].y,
                v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                v_color_unicurvature[v].y,
                v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
                v_color_curvature[v].y,
                v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                v_color_gaussian_curv[v].y,
                v_color_gaussian_curv[v].z;
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
                                  Mesh::Vertex_property<Color> color_prop, Scalar min_value,
                                  Scalar max_value, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    if (min_value == 0.0 && max_value == 0.0) {
        // discard upper and lower bound
        unsigned int n = values.size()-1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        min_value = values[i];
        max_value = values[n-1-i];
    }

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}


}


