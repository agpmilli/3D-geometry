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
      points.conservativeResize(2, points.cols() + 1);
      int n = points.cols() - 1;

      // Compute s using epsilon and n
      double s = epsilon * n;

      // phi is the integral of k(s) on s
      double phi = s; // first curvature k(s) = 1;
      //double phi = (double)1/2 * pow(s,2); // second curvature k(s) = s
      //double phi = (double)1/3 * pow(s,3) - 2.19*s; // third curvature k(s) = s^2 - 2.19

      // Compute next point
      points(0, n) = points(0,n-1) + tangent(0);
      points(1, n) = points(1,n-1) + tangent(1);

      // Compute tangent for next point
      tangent=Vec2(cos(phi), sin(phi));

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
    //points.col (0) = Vec2 (0, 0.0); // second
    // points.col (0) = Vec2 (0, 0.0); // third

    tangent = Vec2 (1.0, 0.0);

    epsilon = 0.05;

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
