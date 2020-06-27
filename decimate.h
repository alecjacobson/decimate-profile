#include <Eigen/Core>
#include <igl/igl_inline.h>

#include "is_edge_manifold.h"
#include "min_heap.h"
#include "collapse_edge.h"
#include "edge_flaps.h"

#include <igl/get_seconds.h>
#include <igl/edge_flaps.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice_mask.h>
#include <igl/slice.h>
#include <igl/remove_unreferenced.h>
#include <igl/collapse_edge.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/is_edge_manifold.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/write_triangle_mesh.h>

#include <iostream>
template <bool old>
bool decimate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const size_t max_m,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G)
{
  ce_t_pop = 0;
  ce_t_cir = 0;
  ce_t_up = 0;
  ce_t_re = 0;
  ce_t_collapse = 0;
  const auto & tictoc = []() -> double
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  tictoc();
  // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV OV;
  DerivedF OF;
  igl::connect_boundary_to_infinity(V,F,OV,OF);
  std::cout<<"connect:      "<<tictoc()<<std::endl;

  Eigen::VectorXi I,J;
  bool ret = false;
  {
    // Decimate 1
    using namespace Eigen;
    using namespace std;
    // Working copies
    Eigen::MatrixXd V = OV;
    Eigen::MatrixXi F = OF;
    VectorXi EMAP;
    MatrixXi E,EF,EI;
    ::edge_flaps(F,E,EMAP,EF,EI);
    std::cout<<"edge_flaps:"<<tictoc()<<std::endl;
    {
      Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> BF;
      Eigen::Array<bool,Eigen::Dynamic,1> BE;
      if(!::is_edge_manifold(F,E.rows(),EMAP,BF,BE))
      {
        return false;
      }
    }
    std::cout<<"is_edge_manif:"<<tictoc()<<std::endl;

    const auto always_try = [](
      const Eigen::MatrixXd &                                         ,/*V*/
      const Eigen::MatrixXi &                                         ,/*F*/
      const Eigen::MatrixXi &                                         ,/*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,/*EF*/
      const Eigen::MatrixXi &                                         ,/*EI*/
      const igl::min_heap< std::tuple<double,int,int> > &             ,/*Q*/
      const Eigen::VectorXi &                                         ,/*EQ*/
      const Eigen::MatrixXd &                                         ,/*C*/
      const int                                                        /*e*/
      ) -> bool { return true;};
    const auto never_care = [](
      const Eigen::MatrixXd &                                         ,   /*V*/
      const Eigen::MatrixXi &                                         ,   /*F*/
      const Eigen::MatrixXi &                                         ,   /*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,  /*EF*/
      const Eigen::MatrixXi &                                         ,  /*EI*/
      const igl::min_heap< std::tuple<double,int,int> > &             ,   /*Q*/
      const Eigen::VectorXi &                                         ,   /*EQ*/
      const Eigen::MatrixXd &                                         ,   /*C*/
      const int                                                       ,   /*e*/
      const int                                                       ,  /*e1*/
      const int                                                       ,  /*e2*/
      const int                                                       ,  /*f1*/
      const int                                                       ,  /*f2*/
      const bool                                                  /*collapsed*/
      )-> void { };
    const auto cost_and_placement = igl::shortest_edge_and_midpoint;
    const auto stopping_condition = [&m,&orig_m,&max_m](
      const Eigen::MatrixXd &                                         ,   /*V*/
      const Eigen::MatrixXi &                                         ,   /*F*/
      const Eigen::MatrixXi &                                         ,   /*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,  /*EF*/
      const Eigen::MatrixXi &                                         ,  /*EI*/
      const igl::min_heap< std::tuple<double,int,int> > &             ,   /*Q*/
      const Eigen::VectorXi &                                         ,   /*EQ*/
      const Eigen::MatrixXd &                                         ,   /*C*/
      const int                                                       ,   /*e*/
      const int                                                       ,  /*e1*/
      const int                                                       ,  /*e2*/
      const int                                                       f1,
      const int                                                       f2 
      )->bool
    {
      // Only subtract if we're collapsing a real face
      if(f1 < orig_m) m-=1;
      if(f2 < orig_m) m-=1;
      return m<=(int)max_m;
    };
    const auto pre_collapse = always_try;
    const auto post_collapse = never_care;

    igl::min_heap<std::tuple<double,int,int> > Q;
    // Could reserve with https://stackoverflow.com/a/29236236/148668
    Eigen::VectorXi EQ = Eigen::VectorXi::Zero(E.rows());
    // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd C(E.rows(),V.cols());
    // Pushing into a vector then using constructor was slower. Maybe using
    // std::move + make_heap would squeeze out something?
    
    tictoc();
    // Separating the cost/placement evaluation from the Q filling is a
    // performance hit for serial but faster if we can parallelize the
    // cost/placement.
    {
      Eigen::VectorXd costs(E.rows());
      igl::parallel_for(E.rows(),[&](const int e)
      {
        double cost = e;
        RowVectorXd p(1,3);
        cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
        C.row(e) = p;
        costs(e) = cost;
      },10000);
      for(int e = 0;e<E.rows();e++)
      {
        Q.emplace(costs(e),e,0);
      }
    }

    std::cout<<"  Qfill: "<<tictoc()<<std::endl;

    int prev_e = -1;
    bool clean_finish = false;

    while(true)
    {
      if(Q.empty())
      {
        break;
      }
      if(std::get<0>(Q.top()) == std::numeric_limits<double>::infinity())
      {
        // min cost edge is infinite cost
        break;
      }
      int e,e1,e2,f1,f2;
      if(collapse_edge<old>(
        cost_and_placement, pre_collapse, post_collapse,
        V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2))
      {
        if(stopping_condition(V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2))
        {
          clean_finish = true;
          break;
        }
      }else
      {
        if(prev_e == e)
        {
          assert(false && "Edge collapse no progress... bad stopping condition?");
          break;
        }
        // Edge was not collapsed... must have been invalid. collapse_edge should
        // have updated its cost to inf... continue
      }
      prev_e = e;
    }
    std::cout<<"    ce_t_pop:"<<ce_t_pop<<std::endl;
    std::cout<<"    ce_t_cir:"<<ce_t_cir<<std::endl;
    std::cout<<"    ce_t_col:"<<ce_t_collapse<<std::endl;
    std::cout<<"    ce_t_up :"<<ce_t_up<<std::endl;
    std::cout<<"    ce_t_re :"<<ce_t_re<<std::endl;
    std::cout<<"  while:"<<tictoc()<<std::endl;
    //igl::write_triangle_mesh("before-ru.obj",V,F);
    // remove all IGL_COLLAPSE_EDGE_NULL faces
    MatrixXi F2(F.rows(),3);
    J.resize(F.rows());
    int m = 0;
    for(int f = 0;f<F.rows();f++)
    {
      if(
        F(f,0) != IGL_COLLAPSE_EDGE_NULL || 
        F(f,1) != IGL_COLLAPSE_EDGE_NULL || 
        F(f,2) != IGL_COLLAPSE_EDGE_NULL)
      {
        F2.row(m) = F.row(f);
        J(m) = f;
        m++;
      }
    }
    std::cout<<"  F2: "<<tictoc()<<std::endl;
    F2.conservativeResize(m,F2.cols());
    J.conservativeResize(m);
    VectorXi _1;
    igl::remove_unreferenced(V,F2,U,G,_1,I);
    ret = clean_finish;
    std::cout<<"  ru: "<<tictoc()<<std::endl;
    //igl::writeDMAT("EQ.dmat",EQ.cast<double>(),false);
  }
  std::cout<<"  clean: "<<tictoc()<<std::endl;
  const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
  igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
  igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
  Eigen::VectorXi _1,I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
  igl::slice(Eigen::VectorXi(I),I2,1,I);
  std::cout<<"postdecimate: "<<tictoc()<<std::endl;
  return ret;
}
