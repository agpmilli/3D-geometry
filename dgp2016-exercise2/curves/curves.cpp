#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
  PointsRenderer render_points = PointsRenderer ();
  SegmentsRenderer render_segments = SegmentsRenderer ();

  MatMxN points;
  MatMxN points_3d_render;
  int num_points;
  SegmentsRenderer::Segments segments;

// ============================================================================
// Exercise 2 : fill the 2 functions below (see PDF for instructions)
// To test your implementation, use the S key for laplacian smoothing and the
// C key for the osculating circle.
// Hint : try to play with epsilon
// ============================================================================
  // time step for smoothing
  double epsilon = 0.001;

  void laplacianSmoothing() {
    // Curve Smoothing - centroid (this function should do one iteration of smoothing)
  }

  void osculatingCircle() {
    // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
  }

// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code)
// ============================================================================

  void generateRandomizedClosedPolyline() {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5*3e-2);

    Vec2 center(3e-2, 2e-3);
    const double radius = 0.3;

    points = MatMxN::Zero(2, num_points);
    for (int i = 0; i < num_points; ++i)
    {
      double frac = static_cast<double>(i) / static_cast<double>(num_points);
      points(0, i) = center(0) + radius * cos (2. * M_PI * frac) + distribution(generator);
      points(1, i) = center(1) + radius * sin (2. * M_PI * frac) + distribution(generator);
    }
  }

  void render () {

    // Prepare the render points
    points_3d_render = MatMxN::Zero(3, points.cols());
    points_3d_render.block(0, 0, 2, points.cols()) = points;

    // Rebuild the segments
    segments.clear();
    for (int i = 0; i < points_3d_render.cols(); ++i) {
      segments.push_back({ points_3d_render.col(i), points_3d_render.col((i+1) % points_3d_render.cols()) });
    }
    render_points.init_data(points_3d_render);
    render_segments.init_data(segments);
  }

  MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
    num_points = 50;
    generateRandomizedClosedPolyline();

    this->scene.add(render_points);
    this->scene.add(render_segments);

    render();
  }

  bool key_callback(int key, int scancode, int action, int mods) override {
    TrackballWindow::key_callback(key, scancode, action, mods);
    if (key == GLFW_KEY_S && action == GLFW_RELEASE)
    {
      laplacianSmoothing();     
    }
    else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
    {
      osculatingCircle();
    }

    render();
    return true;
  }
};


int main(int argc, char** argv)
{
  MainWindow window(argc, argv);
  return window.run();
}
