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
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

#include <cmath>


using std::min;
using std::max;
using namespace surface_mesh;

typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Viewer::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature = mesh.vertex_property<Scalar>("v:unicurvature", 0);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate mean curvature using
    // the length of the uniform Laplacian approximation
    // Save your approximation in unicurvature vertex property of the mesh.

    // iterate over all non-boundary vertices v
    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)){
            unsigned int N = mesh.valence(v);
            Point sum(0.0, 0.0, 0.0);
            // iterate over all neighbors of v
            for(auto neighbor: mesh.vertices(v)){
                sum += (mesh.position(neighbor) - mesh.position(v));
            }
            v_unicurvature[v] = (norm(sum)/(double)N);
        }
    }
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Viewer::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
    Mesh::Vertex_property<Scalar>  v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // run calc_weights function to fill e_weight and v_weight arrays
    calc_weights();
    //iterate over each vertex v
    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)){
            Point sum(0.0, 0.0, 0.0);
            double w = v_weight[v];
            // iterate over each neighbor vi
            for(auto vi: mesh.vertices(v)){
                // get the edge between v and vi and get its weight from e_weight
                auto ei = mesh.find_edge(v, vi);
                double wi = e_weight[ei];
                sum += wi*(mesh.position(vi)-mesh.position(v));
            }
            sum *= w;
            v_curvature[v] = norm(sum);
        }
    }

    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Viewer::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // run calc_weights function to fill e_weight and v_weight arrays
    calc_weights();
    // Iterate over each vertex v
    for(auto v: mesh.vertices()){
        if(!mesh.is_boundary(v)){
            double angle_sum = 0;
            double w = v_weight[v];
            // iterate over all neighbors of v
            for(auto v1: mesh.vertices(v)){
                // get the edge between v and vi and get its weight from e_weight
                auto e1 = mesh.find_edge(v, v1);
                for(auto v2: mesh.vertices(v)){
                    auto e2 = mesh.find_edge(v,v2);
                    if(mesh.find_edge(v1,v2)!=Surface_mesh::Edge()){
                        /* TODO how to get cross product of e1 and e2 ?
                         * angle_sum += acos(e1 * e2);*/
                    }
                }
            }
            v_gauss_curvature[v] = (2*w)*(2.0* M_PI - angle_sum);
        }
    }
    // ------------- IMPLEMENT HERE ---------
}

// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Viewer::uniform_smooth(unsigned int n_iters) {

    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        // ------------- IMPLEMENT HERE ---------
    }

    // update face and vertex normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Viewer::smooth(unsigned int n_iters) {

    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // Perform Laplace-Beltrami smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        // 2) for each non-boundary vertex, update its position using the normalized Laplace-Beltrami operator
        //    (Hint: use the precomputed edge weights in the edge property "e:weight")
        // ------------- IMPLEMENT HERE ---------
    }


    // update face and vertex normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}

// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::uniform_laplacian_enhance_feature(int enhancement_smoothing_iterations,
                                               float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the uniform Laplacian operator:
    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------

    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::laplace_beltrami_enhance_feature(int enhancement_smoothing_iterations,
                                              float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the Laplace-Beltrami operator:
    // 1) perform Laplace-Beltrami smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
