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
    for (int i = 0; i < num_iterations; ++i)
    {
        split_long_edges ();
        collapse_short_edges ();
        equalize_valences ();
        tangential_relaxation ();
    }
}

void MeshProcessing::calc_target_length (const REMESHING_TYPE &remeshing_type) {
    Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator  vv_c, vv_end;
    Scalar                   length;
    Scalar                   mean_length;
    Scalar                   H;
    Scalar                   K;

    Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:meancurvature", 0);
    Mesh::Vertex_property<Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    Mesh::Vertex_property<Scalar> target_new_length  = mesh_.vertex_property<Scalar>("v:newlength", 0);


    if (remeshing_type == AVERAGE)
    {
      // calculate desired length
      for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
        length = 1.0;
        if (!mesh_.is_boundary(*v_it)) {
          //
        }
        target_length[*v_it] = length;
      }

      // smooth desired length
      for (int i = 0; i < 5; i++) {
        //
      }

      // rescale desired length:

    }
    else if (remeshing_type == CURV)
    {

    }
    else if (remeshing_type == HEIGHT)
    {

    }

}

void MeshProcessing::split_long_edges ()
{
    Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
    Mesh::Vertex   v0, v1, v;
    bool            finished;
    int             i;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished=false, i=0; !finished && i<100; ++i)
    {
        //set tu true at the beginning and later set to false if at least one edge is split
        finished = true;
        // Iterate over each edge e
        for(auto e:mesh_.edges()){
            // find the two endpoints of e (v1, v2)
            auto v1 = mesh_.vertex(e, 0);
            auto v2 = mesh_.vertex(e, 1);
            // compute e's edge target length
            Scalar target_length_e = (target_length[v1] + target_length[v2])/2.0;

            //split e if its length is bigger than 4/3 times its target length
            if(mesh_.edge_length(e) > target_length_e*(4.0/3.0)){
                finished = false;
                //add a new vertex v0
                auto v0 = mesh_.add_vertex(Point((mesh_.position(v1)+mesh_.position(v2))/2));
                //split e
                mesh_.split(e, v0);
                //compute and store v0's normal
                normals[v0] = mesh_.compute_vertex_normal(v0);
                //update target edge length
                target_length[v0] = (target_length[v1]+target_length[v2])/2.0;
            }
        }
        // since we created vertices and thus new faces, we should update each face's normal
        mesh_.update_face_normals();
    }
}

void MeshProcessing::collapse_short_edges ()
{
    Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
    Mesh::Vertex   v0, v1;
    Mesh::Halfedge  h01, h10;
    bool            finished, b0, b1;
    int             i;
    bool            hcol01, hcol10;

    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished=false, i=0; !finished && i<100; ++i)
    {
        finished = true;
        // TODO

        for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
        {
            if (!mesh_.is_deleted(*e_it)) // might already be deleted
            {
            }
        }
    }

    mesh_.garbage_collection();

    if (i==100) std::cerr << "collapse break\n";
}

void MeshProcessing::equalize_valences ()
{
    Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
    Mesh::Vertex   v0, v1, v2, v3;
    Mesh::Halfedge   h;
    int             val0, val1, val2, val3;
    int             val_opt0, val_opt1, val_opt2, val_opt3;
    int             ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool            finished;
    int             i;


    // flip all edges
    for (finished=false, i=0; !finished && i<100; ++i)
    {
        finished = true;
        // TODO


        for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
        {
            if (!mesh_.is_boundary(*e_it))
            {
            }
        }
    }

    if (i==100) std::cerr << "flip break\n";
}

void MeshProcessing::tangential_relaxation ()
{
    Mesh::Vertex_iterator     v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    int    valence;
    Point     u, n;
    Point     laplace;

    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");


    // smooth
    for (int iters=0; iters<10; ++iters)
    {
        for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
        {
            if (!mesh_.is_boundary(*v_it))
            {
            }
        }

        for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
            if (!mesh_.is_boundary(*v_it))
                mesh_.position(*v_it) += update[*v_it];
    }
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


