#ifndef IGL_MIN_HEAP_H
#define IGL_MIN_HEAP_H
#include <queue>
#include <vector>
#include <functional>
namespace igl
{
  template<class T> using min_heap = 
    std::priority_queue< T, std::vector<T >, std::greater<T > >;
}
#endif 

