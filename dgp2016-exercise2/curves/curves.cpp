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
  double epsilon1 = 0.2;
  double epsilon2 = 0.001;

  void laplacianSmoothing() {
    // Curve Smoothing - centroid (this function should do one iteration of smoothing)

      // Original length
      double L = 0;
      // After transformation length
      double L1 = 0;

    for (int i = 1; i < num_points-1; i++)
    {
      // get size of each line
      L+= sqrt(pow(points(0,i+1)-points(0,i),2) + pow(points(1,i+1)-points(1,i),2));
      // Do the transformation for current point
      points(0, i) = (1-epsilon1) * points(0, i) + (epsilon1*((points(0, i-1)+points(0, i+1))/2));
      points(1, i) = (1-epsilon1) * points(1, i) + (epsilon1*((points(1, i-1)+points(1, i+1))/2));
      // get size of each line after transformation
      L1+= sqrt(pow(points(0,i+1)-points(0,i),2) + pow(points(1,i+1)-points(1,i),2));
    }

    // Special for last and first points (same as in the for loop)
    L+= sqrt(pow(points(0,1)-points(0,0),2) + pow(points(1,1)-points(1,0),2));
    points(0, 0) = (1-epsilon1) * points(0, 0) + (epsilon1*((points(0, num_points-1)+points(0, 1))/2));
    points(1, 0) = (1-epsilon1) * points(1, 0) + (epsilon1*((points(1, num_points-1)+points(1, 1))/2));
    L1+= sqrt(pow(points(0,1)-points(0,0),2) + pow(points(1,1)-points(1,0),2));

    L+= sqrt(pow(points(0,num_points-1)-points(0,num_points-1),2) + pow(points(1,0)-points(1,0),2));
    points(0, num_points-1) = (1-epsilon1) * points(0, num_points-1) + (epsilon1*((points(0, num_points-2)+points(0, 0))/2));
    points(1, num_points-1) = (1-epsilon1) * points(1, num_points-1) + (epsilon1*((points(1, num_points-2)+points(1, 0))/2));
    L1+= sqrt(pow(points(0,num_points-1)-points(0,num_points-1),2) + pow(points(1,0)-points(1,0),2));

    // uniformly scale the curve back to its original length
    std::cerr << "L = " << L << " and L1 = " << L1;
    points = points * (L/L1);

  }

  void osculatingCircle() {
    // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)

      // Original length
      double L = 0;
      // After transformation length
      double L1 = 0;
      double L_after = 0;

      for (int i = 0; i < num_points; i++)
      {
        // get the coordinates of current point
        Vec2 point = Vec2(points(0,i), points(1,i));
        // create variables for coordinates of 3 points used to find the circumscribed circle
        double point0x, point0y, point1x, point1y, point2x, point2y;
        point1x = points(0,i);
        point1y = points(1,i);
        if(i==num_points-1){
            // get size of each line
            L+= sqrt(pow(points(0,0)-points(0,i),2) + pow(points(1,0)-points(1,i),2));
            point2x = points(0,0);
            point2y = points(1,0);
        } else {
            // get size of each line
            L+= sqrt(pow(points(0,i+1)-points(0,i),2) + pow(points(1,i+1)-points(1,i),2));
            point2x = points(0,i+1);
            point2y = points(1,i+1);
        }
        if(i==0){
            point0x = points(0,num_points-1);
            point0y = points(1,num_points-1);
        } else {
            point0x = points(0,i-1);
            point0y = points(1,i-1);
        }

        // find center of the circumscribed circle
        double D = 2 * ((point0x*(point1y-point2y)) + (point1x*(point2y-point0y)) + (point2x*(point0y-point1y)));
        double Cx=(((pow(point0x,2) + pow(point0y,2))*(point1y - point2y)) + ((pow(point1x,2) + pow(point1y,2))*(point2y - point0y)) + ((pow(point2x,2) + pow(point2y,2))*(point0y - point1y))) / D;
        double Cy=(((pow(point0x,2) + pow(point0y,2))*(point2x - point1x)) + ((pow(point1x,2) + pow(point1y,2))*(point0x - point2x)) + ((pow(point2x,2) + pow(point2y,2))*(point1x - point0x))) / D;
        Vec2 C = Vec2(Cx, Cy);
        Vec2 C_P = C-point;
        double C_P_norm = sqrt(pow(C_P(0),2) + pow(C_P(1),2));
        point = point + (epsilon2 * (C_P / pow(C_P_norm,2)));
        points(0, i) = point(0);
        points(1, i) = point(1);

        if(i==num_points-1){
            // get size of each line after transformation
            L1+= sqrt(pow(points(0,0)-points(0,i),2) + pow(points(1,0)-points(1,i),2));
        } else {
            // get size of each line after transformation
            L1+= sqrt(pow(points(0,i+1)-points(0,i),2) + pow(points(1,i+1)-points(1,i),2));
        }

      }
      // uniformly scale the curve back to its original length
      std::cerr << "L = " << L << " and L1 = " << L1 << " L/L1 = " << L/L1;
      //std::cerr << "pointx = " << points(0,0) << " pointy = " << points(1,0);
      points = points * (L/L1);
      //std::cerr << "pointx = " << points(0,0) << " pointy = " << points(1,0);

      // TODO why doesn't it scale back to initial size??? (TODO PAPER WORK)
      for (int i = 0; i < num_points; i++){
          if(i==num_points-1){
              L_after += sqrt(pow(points(0,0)-points(0,i),2) + pow(points(1,0)-points(1,i),2));
          } else {
              L_after += sqrt(pow(points(0,i+1)-points(0,i),2) + pow(points(1,i+1)-points(1,i),2));
          }
      }
      std::cerr << " L after points * L/L1 = " << L_after;
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
