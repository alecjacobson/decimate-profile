#include <igl/decimate.h>
#include <igl/qslim.h>
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/write_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>


#include <igl/collapse_edge.h> //for IGL_COLLAPSE_EDGE_NULL
#include <igl/circulation.h> //for IGL_COLLAPSE_EDGE_NULL
#include <igl/edge_collapse_is_valid.h> //for IGL_COLLAPSE_EDGE_NULL

#include "min_heap.h"
#include <igl/writeDMAT.h>
#include <utility>

double ce_t_pop = 0;
double ce_t_up = 0;
double ce_t_re = 0;
double ce_t_collapse = 0;
bool collapse_edge(
  const std::function<void(
    const int,
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    double &,
    Eigen::RowVectorXd &)> & cost_and_placement,
  const std::function<bool(
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
    )> & pre_collapse,
  const std::function<void(
    const Eigen::MatrixXd &                                         ,   /*V*/
    const Eigen::MatrixXi &                                         ,   /*F*/
    const Eigen::MatrixXi &                                         ,   /*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,  /*EF*/
    const Eigen::MatrixXi &                                         ,  /*EI*/
    const igl::min_heap< std::tuple<double,int,int> > &             ,   /*Q*/
    const Eigen::VectorXi &                                         ,  /*EQ*/
    const Eigen::MatrixXd &                                         ,   /*C*/
    const int                                                       ,   /*e*/
    const int                                                       ,  /*e1*/
    const int                                                       ,  /*e2*/
    const int                                                       ,  /*f1*/
    const int                                                       ,  /*f2*/
    const bool                                                  /*collapsed*/
    )> & post_collapse,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI,
  igl::min_heap< std::tuple<double,int,int> > & Q,
  Eigen::VectorXi & EQ,
  Eigen::MatrixXd & C,
  int & e,
  int & e1,
  int & e2,
  int & f1,
  int & f2)
{
  const auto & tictoc = []() -> double
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  
  using namespace Eigen;
  using namespace igl;
  std::tuple<double,int,int> p;
  tictoc();
  while(true)
  {
    // Check if Q is empty
    if(Q.empty())
    {
      // no edges to collapse
      return false;
    }
    // pop from Q
    p = Q.top();
    if(std::get<0>(p) == std::numeric_limits<double>::infinity())
    {
      // min cost edge is infinite cost
      return false;
    }
    Q.pop();
    e = std::get<1>(p);
    // Check if matches timestamp
    if(std::get<2>(p) == EQ(e))
    {
      break;
    }
    // must be stale or dead.
    assert(std::get<2>(p)  < EQ(e) || EQ(e) == -1);
    // try again.
  }
  ce_t_pop += tictoc();

  std::vector<int> N  = circulation(e, true,EMAP,EF,EI);
  std::vector<int> Nd = circulation(e,false,EMAP,EF,EI);
  N.insert(N.begin(),Nd.begin(),Nd.end());
  bool collapsed = true;
  if(pre_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e))
  {
    collapsed = collapse_edge(e,C.row(e),V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
  }else
  {
    // Aborted by pre collapse callback
    collapsed = false;
  }
  post_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2,collapsed);
  ce_t_collapse += tictoc();
  if(collapsed)
  {
    // Erase the two, other collapsed edges by marking their timestamps as -1
    EQ(e1) = -1;
    EQ(e2) = -1;
    // update local neighbors
    // loop over original face neighbors
    for(auto n : N)
    {
      if(F(n,0) != IGL_COLLAPSE_EDGE_NULL ||
          F(n,1) != IGL_COLLAPSE_EDGE_NULL ||
          F(n,2) != IGL_COLLAPSE_EDGE_NULL)
      {
        for(int v = 0;v<3;v++)
        {
          // get edge id
          const int ei = EMAP(v*F.rows()+n);
          // compute cost and potential placement
          double cost;
          RowVectorXd place;
          cost_and_placement(ei,V,F,E,EMAP,EF,EI,cost,place);
          // Increment timestamp
          EQ(ei)++;
          // Replace in queue
          Q.emplace(cost,ei,EQ(ei));
          C.row(ei) = place;
        }
      }
    }
    ce_t_up += tictoc();
  }else
  {
    // reinsert with infinite weight (the provided cost function must **not**
    // have given this un-collapsable edge inf cost already)
    // Increment timestamp
    EQ(e)++;
    // Replace in queue
    Q.emplace(std::numeric_limits<double>::infinity(),e,EQ(e));
    ce_t_re += tictoc();
  }
  return collapsed;
}

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
bool decimate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const size_t max_m,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G)
{
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
  if(!igl::is_edge_manifold(OF))
  {
    return false;
  }
  std::cout<<"is_edge_manif:"<<tictoc()<<std::endl;
  Eigen::VectorXi I,J;
  bool ret = false;
  {
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
    // Decimate 1
    using namespace Eigen;
    using namespace std;
    // Working copies
    Eigen::MatrixXd V = OV;
    Eigen::MatrixXi F = OF;
    VectorXi EMAP;
    MatrixXi E,EF,EI;
    igl::edge_flaps(OF,E,EMAP,EF,EI);
    igl::min_heap<std::tuple<double,int,int> > Q;
    // Could reserve with https://stackoverflow.com/a/29236236/148668
    Eigen::VectorXi EQ = Eigen::VectorXi::Zero(E.rows());
    // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd C(E.rows(),V.cols());
    for(int e = 0;e<E.rows();e++)
    {
      double cost = e;
      RowVectorXd p(1,3);
      cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
      C.row(e) = p;
      Q.emplace(cost,e,0);
    }
    int prev_e = -1;
    bool clean_finish = false;
    std::cout<<"  init: "<<tictoc()<<std::endl;

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
      if(collapse_edge(
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
    std::cout<<"    ce_t_col:"<<ce_t_collapse<<std::endl;
    std::cout<<"    ce_t_up :"<<ce_t_up<<std::endl;
    std::cout<<"    ce_t_re :"<<ce_t_re<<std::endl;
    std::cout<<"  while:"<<tictoc()<<std::endl;
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
    F2.conservativeResize(m,F2.cols());
    J.conservativeResize(m);
    VectorXi _1;
    igl::remove_unreferenced(V,F2,U,G,_1,I);
    ret = clean_finish;
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
  ::decimate(V,F,m,U,G);
  std::cout<<"decimate: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("decimate-profile.obj",U,G);
  
  tictoc();
  igl::decimate(V,F,m,U,G,J);
  std::cout<<"igl::decimate: "<<tictoc()<<std::endl;
  igl::write_triangle_mesh("igl-decimate-profile.obj",U,G);
  //tictoc();
  //igl::qslim(V,F,m,U,G,I,J);
  //std::cout<<"   qslim: "<<tictoc()<<std::endl;
  //igl::write_triangle_mesh("qslim-profile.obj",U,G);
}

//#include <queue>
//#include <utility>
//
//int main(int argc, char * argv[])
//{
//  igl::min_heap< std::tuple<double,int,int> > Q;
//}
