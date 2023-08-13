#pragma once


#include "core.hpp"


// #include <chrono>


namespace HW
{

  template<class Node>
  struct default_list_node_traits
  {
    typedef Node node;
    typedef Node* node_ptr;
    typedef const Node* const_node_ptr;
    static node_ptr get_next(const_node_ptr n) { return n->next_; }
    static void set_next(node_ptr n, node_ptr next) { n->next_ = next; }
    static node* get_previous(const_node_ptr n) { return n->prev_; }
    static void set_previous(node_ptr n, node_ptr prev) { n->prev_ = prev; }
  };


  // using namespace std;

  class non_copyable_movable
  {
  protected:
    non_copyable_movable() = default;
    non_copyable_movable(const non_copyable_movable& other) = delete;
    non_copyable_movable(non_copyable_movable&& other) = delete;
    non_copyable_movable& operator=(const non_copyable_movable& other) = delete;
    non_copyable_movable& operator=(non_copyable_movable&& other) = delete;
  };

  

  

  
  




  

  

  

  


  template<class NodeTraits>
  class inorder_iter
  {
    typedef NodeTraits node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    
    node_ptr _node;
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef node_ptr value_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::ptrdiff_t distance_type;
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



  template<typename Algo>
  class tree_inorder_range
  {   
    struct adjusted_traits
    {
      typedef typename Algo::node_ptr node_ptr;

      static node_ptr get_next(node_ptr node) {
        return Algo::next_node(node);
      }
    };


    typedef typename Algo::node_traits node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    
    node_ptr header_;


  public:
    tree_inorder_range(node_ptr header) : header_(header) {}

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
}





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