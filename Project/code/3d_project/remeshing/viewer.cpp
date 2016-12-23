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

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (Screen::keyboardEvent(key, scancode, action, modifiers)) {
        return true;
    }
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        setVisible(false);
        return true;
    }
    return false;
}

void Viewer::draw(NVGcontext *ctx) {
    /* Draw the user interface */
    Screen::draw(ctx);
}

Vector2f Viewer::getScreenCoord() {
    Vector2i pos = mousePos();
    return Vector2f(2.0f * (float)pos.x() / width() - 1.0f,
                    1.0f - 2.0f * (float)pos.y() / height());
}

void Viewer::drawContents() {
    using namespace nanogui;

    /* Draw the window contents using OpenGL */
    shader_.bind();

    Eigen::Matrix4f model, view, proj;
    computeCameraMatrices(model, view, proj);

    Matrix4f mv = view*model;
    Matrix4f p = proj;

    /* MVP uniforms */
    shader_.setUniform("MV", mv);
    shader_.setUniform("P", p);

    // Setup OpenGL (making sure the GUI doesn't disable these
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    // Render everything
    if (wireframe_) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    Vector3f colors(0.98, 0.59, 0.04);
    shader_.setUniform("intensity", colors);
    if (color_mode == CURVATURE) {
        shader_.setUniform("color_mode", int(curvature_type));
    } else {
        shader_.setUniform("color_mode", int(color_mode));
    }
    shader_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());

    if (wireframe_) {
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        colors << 0.0, 0.0, 0.0;
        shader_.setUniform("intensity", colors);
        shader_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    if (normals_) {
        shaderNormals_.bind();
        shaderNormals_.setUniform("MV", mv);
        shaderNormals_.setUniform("P", p);
        shaderNormals_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());
    }
}

bool Viewer::scrollEvent(const Vector2i &p, const Vector2f &rel) {
    if (!Screen::scrollEvent(p, rel)) {
        camera_.zoom = max(0.1, camera_.zoom * (rel.y() > 0 ? 1.1 : 0.9));
    }
    return true;
}

bool Viewer::mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
                              int button, int modifiers) {
    if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
        if (camera_.arcball.motion(p)) {
            //
        } else if (translate_) {
            Eigen::Matrix4f model, view, proj;
            computeCameraMatrices(model, view, proj);
            Point mesh_center = mesh_->get_mesh_center();
            float zval = nanogui::project(Vector3f(mesh_center.x,
                                                   mesh_center.y,
                                                   mesh_center.z),
                                          view * model, proj, mSize).z();
            Eigen::Vector3f pos1 = nanogui::unproject(
                    Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval),
                    view * model, proj, mSize);
            Eigen::Vector3f pos0 = nanogui::unproject(
                    Eigen::Vector3f(translateStart_.x(), mSize.y() -
                                                         translateStart_.y(), zval), view * model, proj, mSize);
            camera_.modelTranslation = camera_.modelTranslation_start + (pos1-pos0);
        }
    }
    return true;
}

bool Viewer::mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
    if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
        if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
            camera_.arcball.button(p, down);
        } else if (button == GLFW_MOUSE_BUTTON_2 ||
                   (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
            camera_.modelTranslation_start = camera_.modelTranslation;
            translate_ = true;
            translateStart_ = p;
        }
    }
    if (button == GLFW_MOUSE_BUTTON_1 && !down) {
        camera_.arcball.button(p, false);
    }
    if (!down) {
        translate_ = false;
    }
    return true;
}

void Viewer::initShaders() {
    // Shaders
    shader_.init(
            "a_simple_shader",

            /* Vertex shader */
            "#version 330\n"
                    "uniform mat4 MV;\n"
                    "uniform mat4 P;\n"
                    "uniform int color_mode;\n"
                    "uniform vec3 intensity;\n"

                    "in vec3 position;\n"
                    "in vec3 valence_color;\n"
                    "in vec3 unicruvature_color;\n"
                    "in vec3 curvature_color;\n"
                    "in vec3 gaussian_curv_color;\n"
                    "in vec3 normal;\n"

                    "out vec3 fcolor;\n"
                    "out vec3 fnormal;\n"
                    "out vec3 view_dir;\n"
                    "out vec3 light_dir;\n"

                    "void main() {\n"
                    "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
                    "    gl_Position = P * vpoint_mv;\n"
                    "    if (color_mode == 1) {\n"
                    "        fcolor = valence_color;\n"
                    "    } else if (color_mode == 2) {\n"
                    "        fcolor = unicruvature_color;\n"
                    "    } else if (color_mode == 3) {\n"
                    "        fcolor = curvature_color;\n"
                    "    } else if (color_mode == 4) {\n"
                    "        fcolor = gaussian_curv_color;\n"
                    "    } else {\n"
                    "        fcolor = intensity;\n"
                    "    }\n"
                    "    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
                    "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
                    "    view_dir = -vpoint_mv.xyz;\n"
                    "}",

            /* Fragment shader */
            "#version 330\n"
                    "uniform int color_mode;\n"
                    "uniform vec3 intensity;\n"

                    "in vec3 fcolor;\n"
                    "in vec3 fnormal;\n"
                    "in vec3 view_dir;\n"
                    "in vec3 light_dir;\n"

                    "out vec4 color;\n"

                    "void main() {\n"
                    "    vec3 c = vec3(0.0);\n"
                    "    if (color_mode == 0) {\n"
                    "        c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
                    "        vec3 n = normalize(fnormal);\n"
                    "        vec3 v = normalize(view_dir);\n"
                    "        vec3 l = normalize(light_dir);\n"
                    "        float lambert = dot(n,l);\n"
                    "        if(lambert > 0.0) {\n"
                    "            c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
                    "            vec3 v = normalize(view_dir);\n"
                    "            vec3 r = reflect(-l,n);\n"
                    "            c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
                    "        }\n"
                    "        c *= fcolor;\n"
                    "    } else {\n"
                    "       c = fcolor;\n"
                    "    }\n"
                    "    if (intensity == vec3(0.0)) {\n"
                    "        c = intensity;\n"
                    "    }\n"
                    "    color = vec4(c, 1.0);\n"
                    "}"
    );

    shaderNormals_.init(
            "normal_shader",
            /* Vertex shader */
            "#version 330\n\n"
                    "in vec3 position;\n"
                    "in vec3 normal;\n"
                    "uniform mat4 MV;\n"
                    "uniform mat4 P;\n"
                    "uniform int normal_selector;\n"
                    "out VS_OUT {\n"
                    "    mat3 normal_mat;\n"
                    "    vec3 normal;\n"
                    "} vs_out;\n"
                    "void main() {\n"
                    "  gl_Position = vec4(position, 1.0);\n"
                    "    vs_out.normal = normal;\n"
                    "    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
                    "}",
            /* Fragment shader */
            "#version 330\n\n"
                    "out vec4 frag_color;\n"
                    "void main() {\n"
                    "   frag_color = vec4(0.0, 1.0, 0.0, 1.0);\n"
                    "}",
            /* Geometry shader */
            "#version 330\n\n"
                    "layout (triangles) in;\n"
                    "layout (line_strip, max_vertices = 6) out;\n"
                    "uniform mat4 MV;\n"
                    "uniform mat4 P;\n"
                    "in VS_OUT {\n"
                    "    mat3 normal_mat;\n"
                    "    vec3 normal;\n"
                    "} gs_in[];\n"
                    "void createline(int index) {\n"
                    "   gl_Position = P * MV * gl_in[index].gl_Position;\n"
                    "   EmitVertex();\n"
                    "   vec4 normal_mv = vec4(normalize(gs_in[index].normal_mat *\n"
                    "                                   gs_in[index].normal), 1.0f);\n"
                    "   gl_Position = P * (MV * gl_in[index].gl_Position\n"
                    "                      + normal_mv * 0.035f);\n"
                    "   EmitVertex();\n"
                    "   EndPrimitive();\n"
                    "}\n"
                    "void main() {\n"
                    "   createline(0);\n"
                    "   createline(1);\n"
                    "   createline(2);\n"
                    "}"
    );
}

Viewer::Viewer() : nanogui::Screen(Eigen::Vector2i(1024, 768), "DGP Viewer") {

    window_ = new Window(this, "Controls");
    window_->setPosition(Vector2i(15, 15));
    window_->setLayout(new GroupLayout());

    PopupButton *popupBtn = new PopupButton(window_, "Open a mesh", ENTYPO_ICON_EXPORT);
    Popup *popup = popupBtn->popup();
    popup->setLayout(new GroupLayout());

    Button* b = new Button(popup, "Bunny");
    b->setCallback([this]() {
        mesh_->load_mesh("../data/bunny.off");
        this->refresh_mesh();
        this->refresh_trackball_center();
    });
    b = new Button(popup, "Max-Planck");
    b->setCallback([this]() {
        mesh_->load_mesh("../data/max.off");
        this->refresh_mesh();
        this->refresh_trackball_center();
    });
    b = new Button(popup, "Geralt");
    b->setCallback([this]() {
        mesh_->load_mesh("../data/geralt.off");
        this->refresh_mesh();
        this->refresh_trackball_center();
    });
    b = new Button(popup, "Geralt no hair");
    b->setCallback([this]() {
        mesh_->load_mesh("../data/geralt_no_hair.off");
        this->refresh_mesh();
        this->refresh_trackball_center();
    });
    b = new Button(popup, "Skull");
    b->setCallback([this]() {
        mesh_->load_mesh("../data/skull.off");
        this->refresh_mesh();
        this->refresh_trackball_center();
    });

    b = new Button(popup, "Open mesh ...");
    b->setCallback([this]() {
        string filename = nanogui::file_dialog({{"obj", "Wavefront OBJ"},
                                         {"ply", "Stanford PLY"},
                                         {"aln", "Aligned point cloud"},
                                         {"off", "Object File Format"}
                                        }, false);
        if (filename != "") {
            mesh_->load_mesh(filename);
            this->refresh_mesh();
            this->refresh_trackball_center();
        }
    });

    b = new Button(popup, "Save current mesh");
    b->setCallback([this]() {
        mesh_->save_mesh();
        this->refresh_mesh();
        this->refresh_trackball_center();
    });

    new Label(window_, "Display Control", "sans-bold");

    b = new Button(window_, "Wireframe");
    b->setFlags(Button::ToggleButton);
    b->setChangeCallback([this](bool wireframe) {
        this->wireframe_ =! this->wireframe_;
    });

    new Label(window_, "Remeshing Type", "sans-bold");
    vector<string> remeshing_type = {"Average", "Curvature", "Height"};
    ComboBox *c = new ComboBox(window_, remeshing_type);
    c->setCallback([this](int remeshing_type) {
        this->remeshing_type = static_cast<mesh_processing::REMESHING_TYPE>
                                                                (remeshing_type);
    });

    b = new Button(window_, "Remeshing");
    b->setCallback([this]() {
        this->mesh_->remesh(this->remeshing_type, 5);
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    new Label(window_, "Smoothing", "sans-bold");
    popupBtn = new PopupButton(window_, "Smooth");
    popup = popupBtn->popup();
    popup->setLayout(new GroupLayout());

    Widget* panelSmooth = new Widget(popup);
    panelSmooth->setLayout(new GroupLayout());

    b = new Button(popup, "Uniform Laplacian");
    b->setCallback([this]() {
        mesh_->uniform_smooth(this->iterationSmoothTextBox->value());
        mesh_->meshProcess();
        this->refresh_mesh();
    });
    b = new Button(popup, "Laplace-Beltrami");
    b->setCallback([this]() {
        mesh_->smooth(this->iterationSmoothTextBox->value());
        mesh_->meshProcess();
        this->refresh_mesh();
    });

    panelSmooth = new Widget(popup);
    GridLayout *layoutSmooth = new GridLayout(Orientation::Horizontal, 2,
                                        Alignment::Middle, 15, 5);
    layoutSmooth->setColAlignment({ Alignment::Maximum, Alignment::Fill });
    layoutSmooth->setSpacing(0, 10);
    panelSmooth->setLayout(layoutSmooth);
    new Label(panelSmooth, "Number of iterations:", "sans-bold");
    iterationSmoothTextBox = new IntBox<int>(panelSmooth, 10);
    iterationSmoothTextBox->setEditable(true);
    iterationSmoothTextBox->setFixedSize(Vector2i(50, 20));
    iterationSmoothTextBox->setDefaultValue("10");
    iterationSmoothTextBox->setFontSize(16);
    iterationSmoothTextBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

    new Label(window_, "Geralt-Effects", "sans-bold");
    popupBtn = new PopupButton(window_, "Mesh Cut");
    popup = popupBtn->popup();
    popup->setLayout(new GroupLayout());

    Widget* panelGeralt = new Widget(popup);
    panelGeralt->setLayout(new GroupLayout());

    b = new Button(panelGeralt, "Separate top of the head (melting-like)");
    b->setCallback([this]() {
        mesh_->separate_head_melting(this->gammaGeraltTextBox->value(), this->lambdaGeraltTextBox->value(), this->thresholdGeraltTextBox->value());
        mesh_->meshProcess();
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    b = new Button(panelGeralt, "Separate top of the head (log curve)");
    b->setCallback([this]() {
        mesh_->separate_head_log(this->gammaGeraltTextBox->value(), this->thresholdGeraltTextBox->value());
        mesh_->meshProcess();
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    b = new Button(panelGeralt, "Create fracture");
    b->setCallback([this]() {
        mesh_->create_fracture();
        mesh_->meshProcess();
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    b = new Button(panelGeralt, "Delete long edges faces");
    b->setCallback([this]() {
        mesh_->delete_long_edges_faces();
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    panelGeralt = new Widget(popup);
    GridLayout *layoutGeralt = new GridLayout(Orientation::Horizontal, 2,
                                        Alignment::Middle, 15, 5);
    layoutGeralt->setColAlignment({ Alignment::Maximum, Alignment::Fill });
    layoutGeralt->setSpacing(0, 10);
    panelGeralt->setLayout(layoutGeralt);
    new Label(panelGeralt, "lambda :", "sans-bold");
    gammaGeraltTextBox = new FloatBox<float>(panelGeralt, 3);
    gammaGeraltTextBox->setEditable(true);
    gammaGeraltTextBox->setFixedSize(Vector2i(50, 20));
    gammaGeraltTextBox->setDefaultValue("3");
    gammaGeraltTextBox->setFontSize(16);
    gammaGeraltTextBox->setFormat("[-]?[0-9]*\\.?[0-9]+");
    new Label(panelGeralt, "threshold :", "sans-bold");
    thresholdGeraltTextBox = new FloatBox<float>(panelGeralt, 20);
    thresholdGeraltTextBox->setEditable(true);
    thresholdGeraltTextBox->setFixedSize(Vector2i(50, 20));
    thresholdGeraltTextBox->setDefaultValue("20");
    thresholdGeraltTextBox->setFontSize(16);
    thresholdGeraltTextBox->setFormat("[-]?[0-9]*\\.?[0-9]+");
    new Label(panelGeralt, "gamma (only used for melting):", "sans-bold");
    lambdaGeraltTextBox = new FloatBox<float>(panelGeralt, 0.2);
    lambdaGeraltTextBox->setEditable(true);
    lambdaGeraltTextBox->setFixedSize(Vector2i(50, 20));
    lambdaGeraltTextBox->setDefaultValue("0.2");
    lambdaGeraltTextBox->setFontSize(16);
    lambdaGeraltTextBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

    new Label(window_, "Skull-Effects", "sans-bold");
    popupBtn = new PopupButton(window_, "Dual Graph");
    popup = popupBtn->popup();
    popup->setLayout(new GroupLayout());

    Widget* panelSkullEffect = new Widget(popup);
    panelSkullEffect->setLayout(new GroupLayout());

    b = new Button(panelSkullEffect, "Make skull pattern with edges");
    b->setCallback([this]() {
        mesh_->skull_dual_graph_pattern();
        mesh_->meshProcess();
        this->mesh_->compute_mesh_properties();
        this->refresh_mesh();
    });

    performLayout();

    initShaders();
    mesh_ = new mesh_processing::MeshProcessing("../data/geralt_no_hair.off");
    this->refresh_mesh();
    this->refresh_trackball_center();
}

void Viewer::refresh_trackball_center() {
    // Re-center the mesh
    Point mesh_center = mesh_->get_mesh_center();
    camera_.arcball = Arcball();
    camera_.arcball.setSize(mSize);
    camera_.modelZoom = 2/mesh_->get_dist_max();
    camera_.modelTranslation = -Vector3f(mesh_center.x, mesh_center.y, mesh_center.z);
}

void Viewer::refresh_mesh() {
    shader_.bind();
    shader_.uploadIndices(*(mesh_->get_indices()));
    shader_.uploadAttrib("position", *(mesh_->get_points()));
    shader_.uploadAttrib("valence_color", *(mesh_->get_colors_valence()));
    shader_.uploadAttrib("unicruvature_color", *(mesh_->get_colors_unicurvature()));
    shader_.uploadAttrib("curvature_color", *(mesh_->get_color_curvature()));
    shader_.uploadAttrib("gaussian_curv_color", *(mesh_->get_colors_gaussian_curv()));
    shader_.uploadAttrib("normal", *(mesh_->get_normals()));
    shader_.setUniform("color_mode", int(color_mode));
    shader_.setUniform("intensity", Vector3f(0.98, 0.59, 0.04));

    shaderNormals_.bind();
    shaderNormals_.shareAttrib(shader_, "indices");
    shaderNormals_.shareAttrib(shader_, "position");
    shaderNormals_.shareAttrib(shader_, "normal");

}

void Viewer::computeCameraMatrices(Eigen::Matrix4f &model,
                                   Eigen::Matrix4f &view,
                                   Eigen::Matrix4f &proj) {

    view = nanogui::lookAt(camera_.eye, camera_.center, camera_.up);

    float fH = std::tan(camera_.viewAngle / 360.0f * M_PI) * camera_.dnear;
    float fW = fH * (float) mSize.x() / (float) mSize.y();

    proj = nanogui::frustum(-fW, fW, -fH, fH, camera_.dnear, camera_.dfar);
    model = camera_.arcball.matrix();

    model = nanogui::scale(model, Eigen::Vector3f::Constant(camera_.zoom * camera_.modelZoom));
    model = nanogui::translate(model, camera_.modelTranslation);
}

Viewer::~Viewer() {
    shader_.free();
    shaderNormals_.free();
}
