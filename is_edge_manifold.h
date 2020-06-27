
#include <igl/unique_edge_map.h>

template <
  typename DerivedF,
  typename DerivedEMAP,
  typename DerivedBF,
  typename DerivedBE>
IGL_INLINE bool is_edge_manifold(
  const Eigen::MatrixBase<DerivedF>& F,
  const typename DerivedF::Index ne,
  const Eigen::MatrixBase<DerivedEMAP>& EMAP,
  Eigen::PlainObjectBase<DerivedBF>& BF,
  Eigen::PlainObjectBase<DerivedBE>& BE)
{
  typedef typename DerivedF::Index Index;
  std::vector<Index> count(ne,0);
  for(Index e = 0;e<EMAP.rows();e++)
  {
    count[EMAP[e]]++;
  }
  const Index m = F.rows();
  BF.resize(m,3);
  BE.resize(ne,1);
  bool all = true;
  for(Index e = 0;e<EMAP.rows();e++)
  {
    const bool manifold = count[EMAP[e]] <= 2;
    all &= BF(e%m,e/m) = manifold;
    BE(EMAP[e]) = manifold;
  }
  return all;
}

template <
  typename DerivedF,
  typename DerivedBF,
  typename DerivedE,
  typename DerivedEMAP,
  typename DerivedBE>
IGL_INLINE bool is_edge_manifold(
  const Eigen::MatrixBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedBF>& BF,
  Eigen::PlainObjectBase<DerivedE>& E,
  Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
  Eigen::PlainObjectBase<DerivedBE>& BE)
{
  using namespace Eigen;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,2> MatrixXF2;
  MatrixXF2 allE;
  igl::unique_edge_map(F,allE,E,EMAP);
  return is_edge_manifold(F,E.rows(),EMAP,BF,BE);
}

template <typename DerivedF>
IGL_INLINE bool is_edge_manifold(
  const Eigen::MatrixBase<DerivedF>& F)
{
  Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> BF;
  Eigen::Array<bool,Eigen::Dynamic,1> BE;
  Eigen::MatrixXi E;
  Eigen::VectorXi EMAP;
  return is_edge_manifold(F,BF,E,EMAP,BE);
}
