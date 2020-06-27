#if true

#include "decimate.h"

#include <igl/decimate.h>
#include <igl/qslim.h>
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/write_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>

#include <igl/writeDMAT.h>
#include <utility>

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

  printf("Hit [Enter] to continue"); getchar(); printf("\n");

  tictoc();
  ::decimate<true>(V,F,m,U,G);
  std::cout<<"deci_old: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("decimate-profile-old.obj",U,G);

  std::this_thread::sleep_for(std::chrono::seconds(3));


  tictoc();
  ::decimate<false>(V,F,m,U,G);
  std::cout<<"decimate: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("decimate-profile.obj",U,G);

  std::this_thread::sleep_for(std::chrono::seconds(3));
  
  tictoc();
  igl::decimate(V,F,m,U,G,J);
  std::cout<<"igl::decimate: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("igl-decimate-profile.obj",U,G);

  std::this_thread::sleep_for(std::chrono::seconds(3));

  igl::qslim(V,F,m,U,G,I,J);
  std::cout<<"   qslim: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("igl-qslim-profile.obj",U,G);
}

//#include <queue>
//#include <utility>
//
//int main(int argc, char * argv[])
//{
//  igl::min_heap< std::tuple<double,int,int> > Q;
//}

#else

#include <igl/is_edge_manifold.h>
#include <igl/matlab_format.h>
#include "circulation.h"
#include "is_edge_manifold.h"
#include "edge_flaps.h"
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
int main(int argc, char * argv[])
{
  //const auto & tictoc = []() -> double
  //{
  //  static double t_start = igl::get_seconds();
  //  double diff = igl::get_seconds()-t_start;
  //  t_start += diff;
  //  return diff;
  //};
  //Eigen::MatrixXd V;
  //Eigen::MatrixXi F;
  //igl::read_triangle_mesh(argv[1],V,F);
  //{
  //  tictoc();
  //  igl::is_edge_manifold(F);
  //  printf("%g secs\n",tictoc());
  //}
  //{
  //  Eigen::VectorXi EMAP;
  //  Eigen::MatrixXi E,EF,EI;
  //  ::edge_flaps(F,E,EMAP,EF,EI);
  //  tictoc();
  //  Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> BF;
  //  Eigen::Array<bool,Eigen::Dynamic,1> BE;
  //  ::is_edge_manifold(F,E.rows(),EMAP,BF,BE);
  //  printf("%g secs\n",tictoc());
  //}

  Eigen::MatrixXi F(4,3);
  F<<
    0,1,3,
    0,2,1,
    0,3,2,
    1,2,3;
  int e = 0;
  //Eigen::MatrixXi F(18,3);
  //F<<
  //  0,4,1,
  //  4,5,1,
  //  1,5,2,
  //  5,6,2,
  //  2,6,3,
  //  6,7,3,
  //  4,8,5,
  //  8,9,5,
  //  5,9,6,
  //  9,10,6,
  //  6,10,7,
  //  10,11,7,
  //  8,12,9,
  //  12,13,9,
  //  9,13,10,
  //  13,14,10,
  //  10,14,11,
  //  14,15,11;
  //int e = 16;
  Eigen::VectorXi EMAP;
  Eigen::MatrixXi E,EF,EI;
  ::edge_flaps(F,E,EMAP,EF,EI);
  //std::cout<<igl::matlab_format(E,"E")<<std::endl;
  printf("E(%d): %d %d\n",e,E(e,0),E(e,1));
  std::vector<int> Nv,Ne,Nf;
  circulation(e,true,F,EMAP,EF,EI,Ne,Nf,Nv);
  printf("Nv: "); for(auto i : Nv) { std::cout<<i<<" "; } printf("\n");
  printf("Ne: "); for(auto i : Ne) { printf("%d,%d ",E(i,0),E(i,1)); } printf("\n");
  printf("Nf: "); for(auto i : Nf) { std::cout<<i<<" "; } printf("\n");

  printf("\n");
  circulation(e,false,F,EMAP,EF,EI,Ne,Nf,Nv);
  printf("Nv: "); for(auto i : Nv) { std::cout<<i<<" "; } printf("\n");
  printf("Ne: "); for(auto i : Ne) { printf("%d,%d ",E(i,0),E(i,1)); } printf("\n");
  printf("Nf: "); for(auto i : Nf) { std::cout<<i<<" "; } printf("\n");

}
#endif
