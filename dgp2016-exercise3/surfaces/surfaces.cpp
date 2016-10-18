//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#include "viewer.h"

using namespace surface_mesh;

// ---------------------------------------- Exercice 2.1 ----------------------------------------

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
        do{
            ++valence;
        }while(++vc != vc_end);
        //store v's valence in the array
        vertex_valence[v] = valence;

    }
}

// ---------------------------------------- Exercice 2.2 ----------------------------------------

// input: vector 3 corners (Vertices) of a triangle
// output: area of the triangle
double computeArea(std::vector<Surface_mesh::Vertex> corners, Surface_mesh *mesh){
    assert(corners.size() == 3);
    // compute the position of the three corners, then the vectors AB and AC of the edges
    auto a = mesh->position(corners[0]);
    auto b = mesh->position(corners[1]);
    auto c = mesh->position(corners[2]);
    auto AB = b-a;
    auto AC = c-a;
    // half cross-product formula: gives the area of a triangle ABC given vectors AB and AC
    double area = 0.5 * sqrt(pow(AB[1]*AC[2]-AB[2]*AC[1], 2) + pow(AB[2]*AC[0]-AB[0]*AC[2],2) + pow(AB[0]*AC[1]-AB[1]*AC[0],2));
    return area;
}

//input: surface mesh, vertex v (the one at the center of the neighborhood), corners contains v and the two other vertices of the triangles
//output: angle at vertex v
double computeAngle(Surface_mesh::Vertex v, std::vector<Surface_mesh::Vertex> corners, Surface_mesh *mesh){
    assert(corners.size() == 3);
    // first we need to figure out which of the three corners is the vertex at the center (call it v_a) of the 1-ring neighborhood:
    unsigned int checkSum = 0;
    Point vertex_A, vertex_B, vertex_C;
    for(int i = 0; i < corners.size(); i++){
        if(corners[i] == v){
            vertex_A = mesh->position(corners[i]);
            vertex_B = mesh->position(corners[(i+1)%3]);
            vertex_C = mesh->position(corners[(i+2)%3]);
            checkSum++;
        }
    }
    assert(checkSum == 1);

    //Then compute the size of the edges a,b,c (a is the edge opposed to v_a, etc)
    auto c = vertex_B - vertex_A; double length_c = sqrt(pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2));
    auto b = vertex_C - vertex_A; double length_b = sqrt(pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2));
    auto a = vertex_C - vertex_B; double length_a = sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));

    // Then compute the angle theta according to the law of cosines:
    double theta = acos((-pow(length_a,2)+pow(length_b,2)+pow(length_c,2))/(2*length_b*length_c));
    return theta;

}

void computeNormalsWithConstantWeights(Surface_mesh *mesh) {
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_cste_weights_n =
            mesh->vertex_property<Point>("v:cste_weights_n", default_normal);
    // initialize constant weight
    unsigned int weight = 1;
    // iterate over all vertices
    for(Surface_mesh::Vertex v: mesh->vertices()){
        // iterate over all triangles incident to v and calculate normal vector multiplied by weight
        Normal sum_n_T(0.0,0.0,0.0);
        for(Surface_mesh::Face f: mesh->faces(v)){
            sum_n_T = sum_n_T + weight * (mesh->compute_face_normal(f));
        }
        auto n_v = sum_n_T.normalize();
        // store vertex normal in the array
        v_cste_weights_n[v] = n_v;
    }
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
    // iterate over all vertices
    for(Surface_mesh::Vertex v: mesh->vertices()){
        // iterate over all triangles incident to v and calculate normal vector multiplied by weight
        Normal sum_n_T((0.0, 0.0, 0.0));
        for(Surface_mesh::Face f: mesh->faces(v)){
            //initialize empty array
            std::vector<Surface_mesh::Vertex> corners;
            // store the three corners of the triangle and use them to compute the area
            for(auto triangle_corner:  mesh->vertices(f)){
                corners.push_back(triangle_corner);
            }
            double weight = computeArea(corners, mesh);
            sum_n_T = sum_n_T + weight * (mesh->compute_face_normal(f));
        }
        auto n_v = sum_n_T.normalize();
        // store vertex normal in the array
        v_area_weights_n[v] = n_v;
    }
}

void computeNormalsWithAngleWeights(Surface_mesh *mesh) {
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_angle_weights_n =
            mesh->vertex_property<Point>("v:angle_weight_n", default_normal);

    // TODO TASK 4
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the angles technique (see .pdf) and store it inside v_angle_weights_n[v]
    for(Surface_mesh::Vertex v: mesh->vertices()){
        // iterate over all triangles incident to v and calculate normal vector multiplied by weight
        Normal sum_n_T(default_normal);
        for(Surface_mesh::Face f: mesh->faces(v)){
            //initialize empty array
            std::vector<Surface_mesh::Vertex> corners;
            // store the three corners of the triangle and use them to compute the area
            for(auto triangle_corner:  mesh->vertices(f)){
                corners.push_back(triangle_corner);
            }
            double weight = computeAngle(v, corners, mesh);
            sum_n_T = sum_n_T + weight * (mesh->compute_face_normal(f));
        }
        auto n_v = sum_n_T.normalize();
        // store vertex normal in the array
        v_angle_weights_n[v] = n_v;

    }
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
