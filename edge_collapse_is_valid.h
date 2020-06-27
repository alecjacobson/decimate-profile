// SIDE EFFECT Nsv,Ndv will be sorted
IGL_INLINE bool edge_collapse_is_valid(
  std::vector<int> & Nsv,
  std::vector<int> & Ndv)
{
  // Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?

  if(Nsv.size()<2 || Ndv.size()<2)
  {
    // Bogus data
    assert(false);
    return false;
  }
  // determine if the first two vertices are the same before reordering.
  // If they are and there are 3 each, then (I claim) this is an edge on a
  // single tet.
  const bool first_two_same = (Nsv[0] == Ndv[0]) && (Nsv[1] == Ndv[1]);
  if(Nsv.size() == 3 && Ndv.size() == 3 && first_two_same)
  {
    // single tet
    return false;
  }
  // https://stackoverflow.com/a/19483741/148668
  std::sort(Nsv.begin(), Nsv.end());
  std::sort(Ndv.begin(), Ndv.end());
  std::vector<int> Nint;
  std::set_intersection(
    Nsv.begin(), Nsv.end(), Ndv.begin(), Ndv.end(), std::back_inserter(Nint));
  // check if edge collapse is valid: intersection of vertex neighbors of s and
  // d should be exactly 2+(s,d) = 4
  // http://stackoverflow.com/a/27049418/148668
  if(Nint.size() != 2)
  {
    return false;
  }
  
  return true;
}
