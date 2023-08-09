#pragma once


#include "core.hpp"


// #include <chrono>


namespace HW
{



//  template<class Real>
//  struct less_real_number		
//  {
//    const Real Epsilon = 1e-20;
//
//    bool operator()(const Real& l, const Real& r) const {
//      return !(abs(r - l) < Epsilon) && l < r;
//    }
//  };

  using namespace std;











//  template<typename Integral>
//  struct integral_forward_iterator
//  {              
//    typedef integral_forward_iterator iterator;
//    typedef Integral integral;
//    typedef Integral* integral_ptr;
//    integral _value;
//
//
//    typedef forward_iterator_tag iterator_category;
//    typedef integral value_type;
//    typedef ptrdiff_t difference_type;
//    typedef ptrdiff_t distance_type;
//    typedef integral_ptr pointer;
//    typedef integral& reference;
//    
//    
//    integral_forward_iterator() {}
//
//    explicit integral_forward_iterator(Integral value) : _value(value) {}
//
//    integral operator*() const {
//      return _value;   
//    }
//
//    integral_ptr operator->() const { 
//      return &_value;
//    }
//    
//    iterator& operator++() {            
//      ++_value;
//      return *this;
//    }
//
//    iterator operator++(int) {
//      integral_forward_iterator result(_value);
//      operator++();
//      return result;
//    }
//
//    friend bool operator==(const iterator& l, const iterator& r) {
//      return l._value == r._value;
//    }
//
//    friend bool operator!=(const iterator& l, const iterator& r) {
//      return !operator==(l, r);
//    }
//  };
//
//
//
//  template<typename Integral>
//  struct range_count_impl
//  {              
//    Integral _count;
//    
//    explicit range_count_impl(const Integral count) : _count(count) {}
//
//    integral_forward_iterator<Integral> begin() {
//      return integral_forward_iterator<Integral>(static_cast<Integral>(0));
//    }
//
//    integral_forward_iterator<Integral> end() {
//      return integral_forward_iterator<Integral>(_count);
//    }    
//  };  
//
//  template<typename Integral>
//  range_count_impl<Integral> range_count(Integral count) {
//    return range_count_impl<Integral>(count);
//  }





//  template 
//  class













  class non_copyable_movable
  {
  protected:
    non_copyable_movable() = default;
    non_copyable_movable(const non_copyable_movable& other) = delete;
    non_copyable_movable(non_copyable_movable&& other) = delete;
    non_copyable_movable& operator=(const non_copyable_movable& other) = delete;
    non_copyable_movable& operator=(non_copyable_movable&& other) = delete;
  };

  template<class T>
  bool has_flag (T value, T flag) {            
    return static_cast<bool>(value & flag);
  }

  template<class T>
  T sqr (T value) {            
    return value*value;
  }

  
  template<class Integral, class Map>
  void for_each(Integral count, Map f){
    for (Integral i = 0; i < count; ++i)
      f(i);
  }

  inline void str_skip_until(char* &ptr, char val) {   
    

    while (*ptr != '\0' && *ptr++ != val) {}    
  }

  inline void str_skip_line(char* &ptr) {   
    str_skip_until(ptr, '\n');
  }

//  inline void str_skip_space(char* &ptr) {       
//    while (*ptr != '\0' && isspace(*ptr))
//      ++ptr;
//  }

  inline void str_skip_word(char* &ptr) {       
    while (isspace(*ptr))
      ++ptr;
    while (*ptr != '\0' && !isspace(*ptr))
      ++ptr;
  }

  template<int Number>
  void str_skip_words(char* &ptr) {       
    for (auto i = 0; i < Number; ++i) {
      while (isspace(*ptr))
        ++ptr;
      while (*ptr != '\0' && !isspace(*ptr))
        ++ptr;
    }
  }

  template<class T>
  void str_parse(char* &ptr, T& val) = delete;

  template<class T>
  T str_parse(char* &ptr) = delete;

  template<>
  inline void str_parse(char* &ptr, unsigned long long& val) {
    val = strtoull(ptr, &ptr, 10);
  }

  template<>
  inline void str_parse(char* &ptr, long& val) {
    val = strtol(ptr, &ptr, 10);
  }

  //HOTFIX
  template<>
  inline void str_parse(char* &ptr, int& val) {
    val = strtol(ptr, &ptr, 10);
  }

  template<>
  inline void str_parse(char* &ptr, unsigned long& val) {
    val = strtoul(ptr, &ptr, 10);
  }

  template<>
  inline void str_parse(char* &ptr, unsigned int& val) {
    val = strtoul(ptr, &ptr, 10);
  }

  template<>
  inline void str_parse(char* &ptr, float& val) {
    val = strtof(ptr, &ptr);    
  }

  template<>
  inline float str_parse(char* &ptr) {
    return strtof(ptr, &ptr);    
  }
  
  template<>
  inline void str_parse(char* &ptr, double& val) {
    val = strtod(ptr, &ptr);    
  }

  template<int Number, class T>
  void str_parse(char* &ptr, T* val) {
    for (auto i = 0; i < Number; ++i)
      str_parse(ptr, *val++);
  }





  

  struct default_avl_balance_node_traits
  {
    typedef dpl::graph::avl_balance balance;

    static balance negative() {
      return dpl::graph::avl_balance::negative_t;
    }

    static balance zero() {
      return dpl::graph::avl_balance::zero_t;
    }

    static balance positive() {
      return dpl::graph::avl_balance::positive_t;
    }
  };

  template<class Node>
  struct default_avl_lrpb_node_traits : default_avl_balance_node_traits // lrpb = left, right, parent, balance
  {  
    using node = Node;
    using node_ptr = node*;
    using const_node_ptr = const node*;
    // TODO!!!!

    static node* get_left(const node* n) {
      return n->left;
    }

    static void set_left(node* n, node* l) {
      n->left = l;
    }

    static node* get_right(const node* n) {
      return n->right;
    }

    static void set_right(node* n, node* r) {
      n->right = r;
    }
    
    static node* get_parent(const node* n) {
      return n->parent;
    }

    static void set_parent(node* n, node* p) {
      n->parent = p;
    }
  };


  template<class NodeTraits>
  class inorder_iter
  {
    typedef NodeTraits node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    
    node_ptr _node;
  public:
    typedef forward_iterator_tag iterator_category;
    typedef node_ptr value_type;
    typedef ptrdiff_t difference_type;
    typedef ptrdiff_t distance_type;
    typedef node_ptr pointer;
    typedef node_ptr reference;
    
    
    inorder_iter() {}

    explicit inorder_iter(node_ptr node)
      : _node(node) {}

    node_ptr operator*() const {
      return _node;   
    }

    node_ptr operator->() const { 
      return _node;
    }


    inorder_iter& operator++() {      
      _node = node_traits::get_next(_node);
      return *this;
    }

    inorder_iter operator++(int) {
      inorder_iter result(_node);
      operator++();
      return result;
    }

    friend bool operator==(const inorder_iter& l, const inorder_iter& r) {
      return l._node == r._node;
    }

    friend bool operator!=(const inorder_iter& l, const inorder_iter& r) {
      return !operator==(l, r);
    }
  };



  template<class NodeAlgorithms>
  class tree_inorder_range
  {   
    struct adjusted_traits
    {
      typedef typename NodeAlgorithms::node_ptr node_ptr;

      static node_ptr get_next(const node_ptr& node) {
        return NodeAlgorithms::next_node(node);
      }
    };


    typedef typename NodeAlgorithms::node_traits node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    
    node_ptr header_;


  public:
    explicit tree_inorder_range(const node_ptr& header) : header_(header) {}

    inorder_iter<adjusted_traits> begin() {
      return inorder_iter<adjusted_traits>(node_traits::get_left(header_));
    }

    inorder_iter<adjusted_traits> end() {
      return inorder_iter<adjusted_traits>(header_);
    }    
  };





 


  
  template <class Integral, class BitType, size_t NumBits>
  struct intergral_plus_shifted_least_significant_bits
  {
    static const auto BitMap = (Integral(1u) << NumBits) - 1;

    static Integral get_value(Integral n) {
      return n >> NumBits;
    }

    static void set_value(Integral& n, Integral p) {
      n = (p << NumBits) | get_bits(n);
    }

    static BitType get_bits(Integral n) {
      return static_cast<BitType>(n & BitMap);
    }

    static void set_bits(Integral& n, BitType c) {
      n = (n & ~BitMap) | static_cast<Integral>(c);
    }

//    static void compress_value(Integral& n) {
//      set_value(n, n);
//    }
//
//    static void decompress_value(Integral& n) {
//      n = get_value(n);
//    }
  };

//  template <class Integral, class BitType, size_t NumBits>
//  struct intergral_plus_aligned_most_significant_bits
//  {
//    
//    static const auto BitShift = (sizeof(Integral)*8 - NumBits);
//    static const auto BitMap = ((size_t(1u) << NumBits) - 1) << BitShift;
//
//    static Integral get_value(Integral n) {
//      return n & ~BitMap;
//    }
//
//    static void set_value(Integral& n, Integral p) {      
//      n = (n & BitMap) | p;
//    }
//
//    static BitType get_shifted_bits(Integral n) {
//      return static_cast<BitType>(n & BitMap);
//    }
//
//    static void set_shifted_bits(Integral& n, BitType c) {      
//      n = (static_cast<Integral>(c)) | (n & ~BitMap);
//    }
//
//
//      //    static BitType get_bits(Integral n) {
//      //     return static_cast<BitType>((n & BitMap) >> BitShift);
//      //   }
//      //
//      //   static void set_bits(Integral& n, BitType c) {      
//      //     n = (static_cast<Integral>(c) << BitShift) | (n & ~BitMap);
//      //   }
//  };


  // NumBits < 3 for x64
//  template <class Pointer, class BitType, size_t NumBits>
//  struct pointer_plus_aligned_least_significant_bits
//  {
//  private:
//    static const auto tag_bitmap_ = (static_cast<size_t>(1u) << NumBits) - 1;
//    static const auto ptr_bitmap_ = ~tag_bitmap_;
//
//  public:
//    typedef Pointer pointer;
//    
//    static pointer get_pointer(pointer n) {
//      return reinterpret_cast<pointer>(reinterpret_cast<size_t>(n) & ptr_bitmap_);
//    }
//
//    static void set_pointer(pointer& n, pointer p) {
//      n = reinterpret_cast<pointer>(reinterpret_cast<size_t>(p) | reinterpret_cast<size_t>(n) & tag_bitmap_);      
//    }
//
//    static BitType get_bits(pointer n) {
//      return static_cast<BitType>(reinterpret_cast<size_t>(n) & tag_bitmap_);      
//    }
//
//    static void set_bits(pointer& n, BitType c) {
//      n = pointer(reinterpret_cast<size_t>(n) & ptr_bitmap_ | static_cast<size_t>(c));
//    }
//
//    static void set_pointer_and_bits(pointer& n, pointer p, BitType c) {
//      n = pointer(reinterpret_cast<size_t>(p) | c);
//    }
//  };

  // NumBits < 3 for x64
  template <class BitType, size_t NumBits>
  struct tagged_pointer_as_size_t
  {
  private:
    static const auto tag_bitmap_ = (static_cast<size_t>(1u) << NumBits) - 1;
    static const auto ptr_bitmap_ = ~tag_bitmap_;

  public:
//    typedef Pointer pointer;
    
    template <class Ptr>
    static Ptr get_pointer(size_t n) {
      return reinterpret_cast<Ptr>(n & ptr_bitmap_);
    }

    template <class Ptr>
    static void set_pointer(size_t& n, Ptr p) {
      n = reinterpret_cast<size_t>(p) | n & tag_bitmap_;      
    }

    static BitType get_bits(size_t n) {
      return static_cast<BitType>(n & tag_bitmap_);      
    }

    static void set_bits(size_t& n, BitType c) {
      n = n & ptr_bitmap_ | static_cast<size_t>(c);
    }

    template <class Ptr>
    static void set_pointer_and_bits(size_t& n, Ptr p, BitType c) {
      n = reinterpret_cast<size_t>(p) | c;
    }
  };



  // 
//  template <class Pointer, class BitType, size_t NumBits>
//  struct pointer_plus_aligned_most_significant_bits
//  {
//    typedef Pointer pointer;
//
//    
//    static const auto BitShift = (sizeof(pointer)*8 - NumBits);
//    static const auto BitMap = ((size_t(1u) << NumBits) - 1) << BitShift;
//
//    static pointer get_pointer(pointer n) {
//      return pointer(size_t(n) & ~BitMap);
//    }
//
//    static void set_pointer(pointer& n, pointer p) {
//      n = pointer(size_t(p) | (size_t(n) & BitMap));      
//    }
//
//    static BitType get_bits(pointer n) {
//      return static_cast<BitType>((size_t(n) & BitMap) >> BitShift);      
//    }
//
//    static void set_bits(pointer& n, BitType c) {      
//      n = pointer(size_t(get_pointer(n)) | size_t(c) << BitShift);
//    }
//
//    static void set_pointer_and_bits(pointer& n, pointer p, BitType c) {
//      n = pointer(size_t(p) | size_t(c) << BitShift);
//    }
//
////    static void set_pointer_and_bits(pointer& n, pointer p, BitType c) {
////      n = pointer(size_t(p) | c);
////    }
//  };











  template <class Container>
  struct container_traits {
    typedef typename Container::value_type value_type;
//    typedef typename Container::pointer pointer;

    static value_type* begin_pointer(Container& c) {
      return &*c.begin();      
    }    
  };

  template <class Array>
  struct container_traits<Array*> {
//    typedef Array value_type;
    static Array* begin_pointer(Array* c) {
      return c;
    }
  };


  template <class Integral>
  Integral integral_power_logN(Integral x, Integral n) {
    if (n == 0)
      return 1;

    Integral y = 1;

    while (n > 1)
      if (n % 2 == 0) {
        x *= x;
        n /= 2;
      }
      else {
        y *= x;
        x *= x;
        n = (n - 1) / 2;
      }

    return x * y;
  }

  template <class Base, class Integral>
  Base integral_power(Base b, Integral n) {
    Base y = 1;

    for (auto i = 0; i < n; i++)
      y *= b;

    return y;
  }

  template <class Integral>
  Integral binary_power(Integral n) {
    return integral_power(static_cast<Integral>(2), n);
  }

  template <class Stream, class Type>
  void write_to_stream(Stream& stream, Type& value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(Type));
  }

  template <class Stream, class Type>
  void read_from_stream(Stream& stream, Type& value) {
    stream.read(reinterpret_cast<char*>(&value), sizeof(Type));
  }


  template <class Stream, class Container>
  void write_to_stream(Stream& stream, Container& container, size_t size) {
    typedef container_traits<Container> traits;            
    stream.write(reinterpret_cast<char*>(traits::begin_pointer(container)), size*sizeof(traits::value_type));
  }

  template <class Stream, class Container>
  void read_from_stream(Stream& stream, Container& container, size_t size) {
    typedef container_traits<Container> traits;            
    container.resize(size);
    stream.read(reinterpret_cast<char*>(traits::begin_pointer(container)), size*sizeof(traits::value_type));
  }

  // inline std::string get_current_executable_path() {
  //   char buffer[MAX_PATH];
  //   GetModuleFileNameA(nullptr, buffer, MAX_PATH);
  //   return std::string(buffer);
  // }

//  template<class T, class I>
//  void reset(unique_ptr<T[]>& uniquePtr, I size) {
//    uniquePtr.reset(new T[size]);
//  }

  


  

  
  template<class T = chrono::milliseconds>
  size_t time_span(chrono::system_clock::time_point t0, chrono::system_clock::time_point t1) {
    return chrono::duration_cast<T>(t1 - t0).count();
  }


  template<class T = chrono::milliseconds>
  size_t time_span(chrono::system_clock::time_point t0) {
    return chrono::duration_cast<T>(chrono::system_clock::now() - t0).count();
  }

//  template<class Container>
//  typename Container::reference push_back_empty(Container& c) {
//    c.push_back(typename Container::value_type());
//    return c.back();
//  }
}
