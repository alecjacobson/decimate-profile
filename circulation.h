//
// for e --> (bf) and ccw=true
//
//     c---d
//    / \ / \
//   a---b-e-f
//    \ / \ /
//     g---h
//
//  Ne = […] -> [fd db dc cb ca ab ag gb gh hb hf fb]
//              {upto cylic order}
//  Nf = […] -> [{bfd}, {bdc}, {bca}, {bag}, {bgh}, {bhf}]
//  Nv = [d c a g h f]
//
IGL_INLINE void circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  std::vector<int> & Ne,
  std::vector<int> & Nf,
  std::vector<int> & Nv)
{
  // prepare output
  Ne.clear();Ne.reserve(10);
  Nv.clear();Nv.reserve(10);
  Nf.clear();Nf.reserve(10);
  const int m = EMAP.size()/3;
  assert(m*3 == EMAP.size());
  const auto & step = [&](
    const int e, 
    const int ff,
    int & ne, 
    int & re,
    int & rv,
    int & nf)
  {
    assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
    //const int fside = EF(e,1)==ff?1:0;
    const int nside = EF(e,0)==ff?1:0;
    const int nv = EI(e,nside);
    // get next face
    nf = EF(e,nside);
    // get next edge 
    const int dir = ccw?-1:1;
    rv = F(nf,nv);
    ne = EMAP(nf+m*((nv+dir+3)%3));
    re = EMAP(nf+m*((nv+2*dir+3)%3));
  };
  // Always start with first face (ccw in step will be sure to turn right
  // direction)
  const int f0 = EF(e,0);
  int fi = f0;
  int ei = e;
  while(true)
  {
    int re,rv;
    step(ei,fi,ei,re,rv,fi);
    Nf.push_back(fi);
    Ne.push_back(re);
    Ne.push_back(ei);
    Nv.push_back(rv);
    // back to start?
    if(fi == f0)
    {
      assert(ei == e);
      break;
    }
  }
}
