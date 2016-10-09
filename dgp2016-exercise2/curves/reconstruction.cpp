#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

# define M_PI           3.14159265358979323846  /* pi */

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
  PointsRenderer render_points = PointsRenderer();
  SegmentsRenderer render_segments = SegmentsRenderer();

  double epsilon;
  MatMxN points;
  MatMxN points_3d_render;
  Vec2 tangent;
  SegmentsRenderer::Segments segments;


// ============================================================================
// Exercise 3 : fill the function below (see PDF for instructions)
// To test your implementation, use the SPACE key
// ============================================================================
  void reconstructCurveStep() {
      //every time you press space, a new point of the curve is added, so the number of points increases
      //std::cerr << "points = " << points << "\n";
      points.conservativeResize(2, points.cols() + 1);
      int n = points.cols() - 1;

      double curvature = 1; // first curve
      //double curvature = epsilon; // first curve
      //double curvature = (pow(s,2)-2.19); // first curve

      // angle = integral of k(s) on s + angle0
      double angle = atan((points(1,n-1) - tangent(1)) /(points(0,n-1) - tangent(0)));
      angle = epsilon + angle*180/M_PI;// first curve k(s) = 1
      // double angle = 1/2 * pow(epsilon,2); // second curve k(s) = s
      // double angle = 1/3 * pow(epsilon,3) - 2,19*epsilon; // third curve k(s) = s^2 - 2.19

      //angle = angle * 180 / M_PI;

      std::cerr << "angle = " << angle << "\n";
      angle = angle*M_PI/180;
      points(0, n) = points(0,n-1) + cos(angle);
      points(1, n) = points(1,n-1) + sin(angle);
      tangent(0) = points(0,n-1) + curvature;
      tangent(1) = points(1,n-1) + curvature;

      std::cerr << "last Point x = " << points(0,n) << "\n";
      std::cerr << "last Point y = " << points(1,n) << "\n";
      std::cerr << "tangent = " << tangent << "\n";
  }
// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code except to uncomment
// the part where you can change the function)
// ============================================================================

  void render() {
    // Prepare the render points
    points_3d_render = MatMxN::Zero(3, points.cols());
    points_3d_render.block(0, 0, 2, points.cols()) = points * 1e-1;

    // Rebuild the segments
    segments.clear();

    for (int i = 1; i < points_3d_render.cols (); ++i) {
      segments.push_back({ points_3d_render.col(i), points_3d_render.col(i-1)});
    }

    render_points.init_data(points_3d_render);
    render_segments.init_data(segments);
  }

  MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
    points = MatMxN(2, 1);

    // First, second, third curve
    points.col(0) = Vec2 (0, -1.0); // first
    // points.col (0) = Vec2 (0, 0.0); // second
    // points.col (0) = Vec2 (0, 0.0); // third

    tangent = Vec2 (1.0, 0.0);

    epsilon = 0.1;

    this->scene.add(render_points);
    this->scene.add(render_segments);

    render ();
  }

  bool key_callback(int key, int scancode, int action, int mods) override {
    TrackballWindow::key_callback(key, scancode, action, mods);
    if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
    {
      reconstructCurveStep();
    }

    render ();
    return true;
  }
};


int main(int argc, char** argv)
{
  MainWindow window(argc, argv);
  return window.run();
}
