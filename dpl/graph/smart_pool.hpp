#pragma once
#include <boost/pool/object_pool.hpp>

namespace HW
{
  // sizeof(T) has to be more than or equal to pointer size (8 bytes in x64)
  template<class T>
  struct smart_pool_traits
  {
     static T* get_next(T* x) {       
       return *reinterpret_cast<T**>(x);
     }

     static void set_next(T* x, T* y) {
       *reinterpret_cast<T**>(x) = y;
     }
  };
  
  template<class T>
  class smart_pool
  { 
    typedef smart_pool_traits<T> traits;

    boost::object_pool<T> _pool;
    T* _releasedStackTop;

  public:
    smart_pool()
      : _releasedStackTop(nullptr) {}

    T* acquire() {
      if (_releasedStackTop) {
        auto result = _releasedStackTop;        
        _releasedStackTop = traits::get_next(_releasedStackTop);       
        return result;
      }      

      return _pool.malloc();
    }

    void release(T* x) {
      traits::set_next(x, _releasedStackTop);
      _releasedStackTop = x;
    }
  };
}
