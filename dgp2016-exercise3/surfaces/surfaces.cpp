//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#include "viewer.h"

using namespace surface_mesh;

void computeValence(Surface_mesh *mesh) {
    Surface_mesh::Vertex_property<unsigned int> vertex_valence =
            mesh->vertex_property<unsigned int>("v:valence", 0);

    // TODO TASK 1
    // compute the valence for each vertex v and store it inside vertex_valence[v]

    //setup vertex iterator
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh->vertices_begin();
    v_end = mesh->vertices_end();

    //setup neighbors iterator
    Surface_mesh::Vertex_around_vertex_circulator vc, vc_end;

    //iterate over all vertices
    for(v_it=v_begin; v_it != v_end; ++v_it){
        Surface_mesh::Vertex v = *v_it;
        //iterate over v's neighbors
        vc = mesh->vertices(v);
        vc_end = vc;
        unsigned int valence = 0;
        //TODO while instead of do while? check
        do{
            ++valence;
        }while(++vc != vc_end);
        //store v's valence in the array
        vertex_valence[v] = valence;

    }


    //for vertex in vertices:
    //valence = 0
    //count the number of neighbors of vertex and assign it to valence
    //vertex_valence[vertex] = valence
}

void computeNormalsWithConstantWeights(Surface_mesh *mesh) {
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_cste_weights_n =
            mesh->vertex_property<Point>("v:cste_weights_n", default_normal);

    // TODO TASK 2
    // Compute the normals for each vertex v in the mesh using the constant weights
    // technique (see .pdf) and store it inside v_cste_weights_n[v]
}

void computeNormalsByAreaWeights(Surface_mesh *mesh) {
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_area_weights_n =
            mesh->vertex_property<Point>("v:area_weight_n", default_normal);

    // TODO TASK 3
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the areas technique (see .pdf) and store inside v_area_weights_n[v]
}

void computeNormalsWithAngleWeights(Surface_mesh *mesh) {
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_angle_weights_n =
            mesh->vertex_property<Point>("v:angle_weight_n", default_normal);

    // TODO TASK 4
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the angles technique (see .pdf) and store it inside v_angle_weights_n[v]
}

// #############################################################################
int main(int /* argc */, char ** /* argv */) {
    try {
        nanogui::init();
        {
            // Load the Mesh
            Surface_mesh mesh;
            if (!mesh.read("../data/bunny.off")) {
                std::cerr << "Mesh not found, exiting." << std::endl;
                return -1;
            }

            computeValence(&mesh);
            computeNormalsWithConstantWeights(&mesh);
            computeNormalsByAreaWeights(&mesh);
            computeNormalsWithAngleWeights(&mesh);

            nanogui::ref<Viewer> app = new Viewer(&mesh);
            app->drawAll();
            app->setVisible(true);
            nanogui::mainloop();
        }

        nanogui::shutdown();
    } catch (const std::runtime_error &e) {
        std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(_WIN32)
        MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
        std::cerr << error_msg << endl;
#endif
        return -1;
    }

    return 0;
}
