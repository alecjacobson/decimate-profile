#include <igl/edge_flaps.h>

IGL_INLINE void edge_flaps(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & uE,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  Eigen::MatrixXi allE;
  igl::unique_edge_map(F,allE,uE,EMAP);
  // Const-ify to call overload
  const auto & cuE = uE;
  const auto & cEMAP = EMAP;
  return igl::edge_flaps(F,cuE,cEMAP,EF,EI);
}

