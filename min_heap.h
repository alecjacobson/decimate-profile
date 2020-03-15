#ifndef IGL_MIN_HEAP_H
#define IGL_MIN_HEAP_H
#include <queue>
#include <vector>
#include <functional>
namespace igl
{
  template<class T> using min_heap = 
    std::priority_queue< T, std::vector<T >, std::greater<T > >;

  // Slow.
  //template <class T>
  //class min_heap : public std::multimap<
  //  typename std::tuple_element<0,T>::type,
  //  std::pair<
  //    typename std::tuple_element<1,T>::type,
  //    typename std::tuple_element<2,T>::type 
  //    >
  //  >
  //{
  //  public:
  //    T top() const 
  //    {
  //      const auto & b = this->begin();
  //      return {
  //        b->first,
  //        b->second.first,
  //        b->second.second
  //      };
  //    };
  //    void pop()
  //    {
  //      this->erase(this->begin());
  //    };
  //    void emplace(
  //      const typename std::tuple_element<0,T>::type & k,
  //      const typename std::tuple_element<1,T>::type & p1,
  //      const typename std::tuple_element<2,T>::type & p2)
  //    {
  //      this->insert({k,{p1,p2}});
  //    };
  //};

  //template <class T>
  //class min_mutable_heap : public std::vector<T>
  //{
  //  public:
  //  const T & top() const
  //  {
  //    return Q.front();
  //  };
  //}
  //void push(const T & val)
  //{
  //  this->push_back(val);
  //  std::push_heap(this->begin(),this->end());
  //};
  //void pop()
  //{
  //  std::pop_heap(this->begin(),this->end());
  //  this->pop_back();
  //};
}
#endif 

