#include <igl/decimate.h>
#include <igl/qslim.h>
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/write_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>

int main(int argc, char * argv[])
{
  const auto & tictoc = []() -> double
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };


  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argv[1],V,F);
  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi I,J;

  const int m = F.rows()*0.1;
  
  tictoc();
  igl::decimate(V,F,m,U,G,J);
  std::cout<<"decimate: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("decimate-profile.obj",U,G);
  tictoc();
  igl::qslim(V,F,m,U,G,I,J);
  std::cout<<"   qslim: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("qslim-profile.obj",U,G);
}
