//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
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

// ============================================================================
// EXERCISE 5.1
// ============================================================================
void MeshProcessing::implicit_smoothing(const double timestep) {

    const int n = mesh_.n_vertices();
    //std::cout << n << std::endl;

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd B(n,3);

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.1 HERE
    // ========================================================================

    for(auto vi: mesh_.vertices()){
        // Get the index and the position of v1
        auto i = vi.idx();
        auto vi_pos = mesh_.position(vi);
        // If v1 is boundary
        if(mesh_.is_boundary(vi)){
            // we put a line in B with the position of vi
            B.row(i)=Eigen::RowVector3d(vi_pos[0], vi_pos[1], vi_pos[2]);
            // and put 1 to its diagonal position in A
            triplets.push_back(Eigen::Triplet<double>(i,i,1));
        }
        // In the other case
        else {
            // We put a line in B with the position of vi multiplied by 2*Ai (multiplication by D^(-1))
            B.row(i)=Eigen::RowVector3d(vi_pos[0], vi_pos[1], vi_pos[2]) / area_inv[vi];
            // and put 2*Ai to its diagonal position in A (multiplication D^(-1))
            triplets.push_back(Eigen::Triplet<double>(i,i, 1/area_inv[vi]));
            // Then we iterate on the neighbors (vjs) of vi
            for(auto vj : mesh_.vertices(vi)){
                // Get the index of current vj
                auto j = vj.idx();
                // We compute Mij of vi and the current vj using timestep and cotan
                auto Mij = timestep * cotan[mesh_.find_edge(vi,vj)];

                // If vj is boundary
                if(mesh_.is_boundary(vj)){
                    auto vj_pos = mesh_.position(vj);
                    // we add to the line already create above the position of vj multiplied by Mij
                    B.row(i) += Eigen::RowVector3d(vj_pos[0], vj_pos[1], vj_pos[2]) * Mij;
                }
                // In the other case
                else {
                    // put -Mij to the position of the intersection between vi and vj (cot(alpha(i,j)) + cot(beta(i,j))
                    triplets.push_back(Eigen::Triplet<double>(i,j, -Mij));
                }
                // We sum the cotan of the neighbors at position (i,i) as seen in the formula (to get - sum(cot(alpha(i,j) + cot(beta(i,j)) with j being the indices of the neighbors of vi (index i))
                triplets.push_back(Eigen::Triplet<double>(i,i,Mij));
            }
        }
    }

    // build sparse matrix from triplets
    A.setFromTriplets(triplets.begin(), triplets.end());

    // solve A*X = B
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(A);
    Eigen::MatrixXd X = solver.solve(B);

    // copy solution
    for (int i = 0; i < n; ++i)
    {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = X(i, dim);
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

// ============================================================================
// EXERCISE 5.2
// ============================================================================
void MeshProcessing::minimal_surface() {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");
    auto points_init = mesh_init_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> L (n, n);
    Eigen::MatrixXd rhs (Eigen::MatrixXd::Zero (n, 3));

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets_L;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.2 HERE
    // ========================================================================
    for(auto vi: mesh_.vertices()){
        // Get the index and the position of vi
        auto i = vi.idx();
        auto vi_pos = mesh_.position(vi);
        // If vi is boundary
        if(mesh_.is_boundary(vi)){
            // we put a line in B with the position of vi
            rhs.row(i)=Eigen::RowVector3d(vi_pos[0], vi_pos[1], vi_pos[2]);
            // and put 1 to its diagonal position in A
            triplets_L.push_back(Eigen::Triplet<double>(i,i,1));
        }
        // In the other case
        else {
            // We iterate on the neighbors (vjs) of vi
            for(auto vj : mesh_.vertices(vi)){
                // Get the index of current vj
                auto j = vj.idx();
                // We compute Mij of vi and the current vj using timestep and cotan
                auto Mij = cotan[mesh_.find_edge(vi,vj)] * area_inv[vi];
                // put -Mij to the position of the intersection between vi and vj (cot(alpha(i,j)) + cot(beta(i,j))
                triplets_L.push_back(Eigen::Triplet<double>(i,j, -Mij));
                // We sum the cotan of the neighbors at position (i,i) as seen in the formula (to get - sum(cot(alpha(i,j) + cot(beta(i,j)) with j being the indices of the neighbors of vi (index i))
                triplets_L.push_back(Eigen::Triplet<double>(i,i,Mij));
            }
        }
    }

    L.setFromTriplets (triplets_L.begin (), triplets_L.end ());

    // solve A*X = B
    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver(L);
    if (solver.info () != Eigen::Success) {
        printf("linear solver init failed.\n");
    }

    Eigen::MatrixXd X = solver.solve(rhs);
    if (solver.info () != Eigen::Success) {
        printf("linear solver failed.\n");
    }

    // copy solution
    for (int i = 0; i < n; ++i) {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim) {
            points[v][dim] += 1. * (X(i, dim) - points[v][dim]);
        }
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}
void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // For each non-boundary vertex, approximate mean curvature using
    // the length of the uniform Laplacian approximation
    // Save your approximation in unicurvature vertex property of the mesh.

    // iterate over all non-boundary vertices v
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

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
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

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
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

void MeshProcessing::uniform_smooth(const unsigned int iterations) {

    for (unsigned int iter=0; iter<iterations; ++iter) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        for(auto v: mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                unsigned int N = mesh_.valence(v);
                Point sum(0.0, 0.0, 0.0);
                // iterate over all neighbors of v
                for(auto neighbor: mesh_.vertices(v)){
                    sum += mesh_.position(neighbor);
                }
                auto vec = (sum/(double)N) - mesh_.position(v);
                // compute the new position of the current vertex
                mesh_.position(v) = mesh_.position(v) + 0.5 * vec;
            }
        }
    }
    // update face and vertex normals
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

void MeshProcessing::smooth(const unsigned int iterations) {
    Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0);
    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
        // Perform Laplace-Beltrami smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        // 2) for each non-boundary vertex, update its position using the normalized Laplace-Beltrami operator
        //    (Hint: use the precomputed edge weights in the edge property "e:weight")

        // precompute edge weights using calc_edge_weights()
        calc_edges_weights();
        for(auto v: mesh_.vertices()){
            if(!mesh_.is_boundary(v)){
                Point sum(0.0, 0.0, 0.0);
                double sum_wi = 0.0;
                // iterate over all neighbors of v
                for(auto neighbor: mesh_.vertices(v)){
                    // get the edge between v and vi and get its weight from e_weight
                    auto ei = mesh_.find_edge(v, neighbor);
                    double wi = e_weight[ei];
                    // compute the sum of wi and the sum wi * (pos(neighbor) - pos(v))
                    sum_wi += wi;
                    sum += wi*(mesh_.position(neighbor)-mesh_.position(v));
                }
                auto vec = 1/sum_wi * sum;
                // compute the new position of the current vertex
                mesh_.position(v) = mesh_.position(v) + 0.5 * vec;
            }
        }
    }
    // update face and vertex normals
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

void MeshProcessing::uniform_laplacian_enhance_feature(const unsigned int iterations,
                                                       const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // Feature enhancement using the uniform Laplacian operator:
    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    std::vector<Point> vs_before;
    for(auto v: mesh_.vertices()){
            // save the position of the vertices before the smoothing
            vs_before.push_back(mesh_.position(v));
    }
    // perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    uniform_smooth(iterations);
    int i = 0;
    for(auto v: mesh_.vertices()){
        // compute the new position of the current vertex
        mesh_.position(v) = mesh_.position(v) + coefficient * (vs_before[i] - mesh_.position(v));
        i+=1;
    }
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

void MeshProcessing::laplace_beltrami_enhance_feature(const unsigned int iterations,
                                                      const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // Feature enhancement using the Laplace-Beltrami operator:
    // 1) perform Laplace-Beltrami smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    std::vector<Point> vs_before;
    for(auto v: mesh_.vertices()){
            // save the position of the vertices before the smoothing
            vs_before.push_back(mesh_.position(v));
    }
    // perform Laplace-Beltrami smoothing for enhancement_smoothing_iterations iterations
    smooth(iterations);
    int i = 0;
    for(auto v: mesh_.vertices()){
            // compute the new position of the current vertex
            mesh_.position(v) = mesh_.position(v) + coefficient * (vs_before[i] - mesh_.position(v));
            i+=1;
    }
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
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
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
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
                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

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

MeshProcessing::~MeshProcessing() {}
}
