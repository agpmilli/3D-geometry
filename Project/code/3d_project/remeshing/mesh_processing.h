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
#ifndef MESH_PROCESSING_H
#define MESH_PROCESSING_H

#include <surface_mesh/Surface_mesh.h>
#include <Eigen/Sparse>

typedef surface_mesh::Surface_mesh Mesh;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;

namespace mesh_processing {

using std::string;
enum REMESHING_TYPE : int { AVERAGE = 0, CURV = 1, HEIGHT = 2 };

class MeshProcessing {

public:
    MeshProcessing(const string& filename);
    ~MeshProcessing();

    const surface_mesh::Point get_mesh_center() { return mesh_center_; }
    const float get_dist_max() { return dist_max_; }
    const Eigen::MatrixXf* get_points() { return &points_; }
    const MatrixXu* get_indices() { return &indices_; }
    const Eigen::MatrixXf* get_normals() { return &normals_; }
    const Eigen::MatrixXf* get_colors_valence() { return &color_valence_; }
    const Eigen::MatrixXf* get_colors_unicurvature() { return &color_unicurvature_; }
    const Eigen::MatrixXf* get_colors_gaussian_curv() { return &color_gaussian_curv_; }
    const Eigen::MatrixXf* get_color_curvature() { return &color_curvature_; }
    const unsigned int get_number_of_face() { return mesh_.n_faces(); }


    void remesh (const REMESHING_TYPE &remeshing_type, const int &num_iterations);
    void calc_target_length (const REMESHING_TYPE &remeshing_type);

    void separate_head();
    void delete_long_edges_faces();
    void delete_big_faces();
    double compute_area_face (Mesh::Face face);

    void make_skull_pattern_edges();
    void create_single_cylinder();
    surface_mesh::Point get_point_from_tuple_vector(Mesh::Face f, std::vector<std::tuple<surface_mesh::Surface_mesh::Face, surface_mesh::Point> > v);
    void MeshProcessing::build_cylinder(surface_mesh::Point p_a, surface_mesh::Point p_b, double r);

    void create_isocahedron(double r, surface_mesh::Point centerPoint);
    surface_mesh::Point middle_point(surface_mesh::Point a, surface_mesh::Point b);
    surface_mesh::Point push_to_radius(surface_mesh::Point point,double radius);
    void create_spheres_on_vertices(std::vector<surface_mesh::Point> dual_intersections);

    void split_long_edges ();
    void collapse_short_edges ();
    void equalize_valences ();
    void tangential_relaxation ();

    void load_mesh(const string& filename);
    void save_mesh();
    void compute_mesh_properties();
    void meshProcess();

    void calc_mean_curvature();
    void calc_uniform_mean_curvature();
    void calc_gauss_curvature();

    void uniform_smooth(unsigned int n_iters);
    void smooth(unsigned int n_iters);

    void uniform_laplacian_enhance_feature(int enhancement_smoothing_iterations, float enhancement_coef);
    void laplace_beltrami_enhance_feature(int enhancement_smoothing_iterations, float enhancement_coef);

private:
    void calc_weights();
    void calc_edges_weights();
    void calc_vertices_weights();

private:
    Mesh mesh_;
    Mesh mesh_init_;
    surface_mesh::Point mesh_center_ = surface_mesh::Point(0.0f, 0.0f, 0.0f);
    float dist_max_ = 0.0f;

    Eigen::MatrixXf points_;
    MatrixXu indices_;
    Eigen::MatrixXf normals_;
    Eigen::MatrixXf color_valence_;
    Eigen::MatrixXf color_unicurvature_;
    Eigen::MatrixXf color_gaussian_curv_;
    Eigen::MatrixXf color_curvature_;

    void color_coding(Mesh::Vertex_property<surface_mesh::Scalar> prop,
                      Mesh *mesh,
                      Mesh::Vertex_property<surface_mesh::Color> color_prop,
                      surface_mesh::Scalar min_value = 0.0,
                      surface_mesh::Scalar max_value = 0.0, int bound = 20);
    void set_color(Mesh::Vertex v, const surface_mesh::Color& col,
                   Mesh::Vertex_property<surface_mesh::Color> color_prop);
    surface_mesh::Color value_to_color(surface_mesh::Scalar value,
                                       surface_mesh::Scalar min_value,
                                       surface_mesh::Scalar max_value);
    };

}

#endif // MESH_PROCESSING_H
