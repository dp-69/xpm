#pragma once

namespace HW { namespace dynamic_connectivity
{ 
  template<class Algorithms>
  class cyclic_operations
  {
    typedef Algorithms algo;    
    typedef typename algo::node_traits node_traits;

    using node = typename node_traits::node;
    using node_ptr = node*;
    using const_node_ptr = const node*;
    
    typedef boost::intrusive::default_path_buffer<node_traits> default_path;

  public:    
    static void principal_cut(node_ptr header, node_ptr least) {      
      if (node_traits::get_left(header) != least) {
        if (node_traits::get_right(header) == least) {
          algo::erase(header, least);
          algo::push_front(header, least);
        }
        else {       
          node header_b_storage;
          auto header_b = &header_b_storage;
          algo::init_header(header_b);          
          algo::split_tree(header, least, header_b);
          algo::push_front(header_b, least);

          auto right_b = node_traits::get_right(header_b);
          algo::erase(header_b, right_b);
          algo::join_trees(header_b, right_b, header);
          algo::swap_tree(header, header_b);
        }
      }
    }

    static void principal_cut_least_dropped(node_ptr header, node_ptr least) {
      if (node_traits::get_right(header) == least || node_traits::get_left(header) == least)
        algo::erase(header, least);
      else {        
        node headerNodeB;
        auto headerB = &headerNodeB;
        algo::init_header(headerB);        

        algo::split_tree(header, least, headerB);

        auto rightB = node_traits::get_right(headerB);
        algo::erase(headerB, rightB);

        algo::join_trees(headerB, rightB, header);

        algo::swap_tree(header, headerB);
      }  
    }

    static node_ptr next_node(node_ptr node) {
      auto next = algo::next_node(node);      
      return node_traits::is_header(next) ? node_traits::get_left(next) : next;       
    }


    // Split initial BST into the following format: A + {n0} + B + {n1} + C.
    // Returns:
    //   headerB <- B (from n0 to n1)
    //   headerA <- C + A (from n1 to n0)

    static void split(node_ptr headerA, node_ptr headerB, node_ptr n0, node_ptr n1) {            
      algo::split_tree(headerA, n0, headerB);
            
      node headerNodeC;
      auto headerC = &headerNodeC;
      algo::init_header(headerC);      

      algo::split_tree(headerB, n1, headerC);
      
      if (node_traits::get_parent(headerA)) {
        auto leftA = node_traits::get_left(headerA);
        algo::erase(headerA, leftA);
        algo::join_trees(headerC, leftA, headerA);
      }
      
      algo::swap_tree(headerA, headerC);      


//      algo::split_tree(headerA, n0, headerB);
//            
//      node headerNodeC;
//      auto headerC = &headerNodeC;
//      algo::init_header(headerC);      
//
//      algo::split_tree(headerB, n1, headerC);
//      
//      algo::swap_tree(headerA, headerC);
//
//      if (node_traits::get_parent(headerC)) {
//        auto leftC = node_traits::get_left(headerC);
//        algo::erase(headerC, leftC);
//        algo::join_trees(headerA, leftC, headerC);
//      }
    }
    


// @@@@ Cut at the first out of the three operands and perform binary comparison on the remainder two.

//    static bool less_than_strict(node_ptr* iter0, node_ptr* iter1, node_ptr* iter2, node_ptr* path0, node_ptr* path1/*, node_ptr* path2*/) {                   
//      auto xSide = algorithms::less_than(iter0, iter1, path0/*, path1*/);
//      auto ySide = algorithms::less_than(iter0, iter2, path0/*, path2*/);
//      return xSide == ySide ? algorithms::less_than(iter1, iter2, path1/*, path2*/) : xSide;            
//    }
//
//    static bool less_than_low1_low2(node_ptr* iter0, node_ptr* iter1, node_ptr* iter2, node_ptr* path0, node_ptr* path1, node_ptr* path2) {                   
//      auto xSide = *(path0 + 1) == *(path1 + 1) || algorithms::less_than(iter0, iter1, path0);
//      auto ySide = *(path0 + 1) == *(path2 + 1) || algorithms::less_than(iter0, iter2, path0);
//      return xSide == ySide ? algorithms::less_than(iter1, iter2, path1) : xSide;            
//    }
//
//    static bool less_than_low1_high2(node_ptr* iter0, node_ptr* iter1, node_ptr* iter2, node_ptr* path0, node_ptr* path1, node_ptr* path2) {                   
//      auto xSide = *(path0 + 1) == *(path1 + 1) || algorithms::less_than(iter0, iter1, path0);
//      auto ySide = *(path0 + 1) == *(path2 + 1) ? false : algorithms::less_than(iter0, iter2, path0);
//      return xSide == ySide ? algorithms::less_than(iter1, iter2, path1) : xSide;            
//    }


    static bool less_than_low_low(const node_ptr& x0, const node_ptr& x1, const node_ptr& x2) {
      auto iter0 = algo::get_path(x0, default_path::path0);
      auto iter1 = algo::get_path(x1, default_path::path1);
      auto iter2 = algo::get_path(x2, default_path::path2);

      auto x1Side = x0 == x1 || algo::less_than(iter0, iter1, default_path::path0);
      auto x2Side = x0 == x2 || algo::less_than(iter0, iter2, default_path::path0);
      return x1Side == x2Side ? algo::less_than(iter1, iter2, default_path::path1) : x1Side;    
    }

//    template<class Tree>
//    class cyclic_comparator_dynamic0_low1_high2
//    {    
//      typedef Tree tree;
//      typedef typename tree::node_traits node_traits;
//      typedef typename node_traits::node_ptr node_ptr;    
//      typedef typename tree::node_algorithms algo;
//      typedef boost::intrusive::default_path_buffer<node_traits> default_path_buffer;
//
//      node_ptr header_;
//      node_ptr x0_;
//      node_ptr* path0_;
//      node_ptr* path1_;
//      node_ptr* path2_;
//
//    public:    
//      cyclic_comparator_dynamic0_low1_high2(const node_ptr& header, const node_ptr& x0, node_ptr* path0, node_ptr* path1, node_ptr* path2)
//        : header_(header), x0_(x0), path0_(path0), path1_(path1), path2_(path2) {}
//
//      explicit cyclic_comparator_dynamic0_low1_high2(const node_ptr& header, const node_ptr& x0)
//        : cyclic_comparator_dynamic0_low1_high2(header, x0, default_path_buffer::path0, default_path_buffer::path1, default_path_buffer::path2) { }
//    
//      bool operator()(const node_ptr& x1, const node_ptr& x2) {                  
//        auto iter0 = algo::get_path(x0_, path0_);
//        auto iter1 = algo::get_path(x1, path1_);
//        auto iter2 = algo::get_path(x2, path2_);
//      
//      
//        auto x1Side = x0_ == x1 || algo::less_than(iter0, iter1, path0_);
//        auto x2Side = x0_ == x2 ? false : algo::less_than(iter0, iter2, path0_);
//        return x1Side == x2Side ? algo::less_than(iter1, iter2, path1_) : x1Side;                              
//      }
//    };





  };

//  template<class Tree>
//  class cyclic_comparator_dynamic0_low1_high2
//  {    
//    typedef Tree tree;
//    typedef typename tree::node_traits node_traits;
//    typedef typename node_traits::node_ptr node_ptr;    
//    typedef typename tree::node_algorithms algo;
//    typedef boost::intrusive::default_path_buffer<node_traits> default_path_buffer;
//
//    node_ptr header_;
//    node_ptr x0_;
//    node_ptr* path0_;
//    node_ptr* path1_;
//    node_ptr* path2_;
//
//  public:    
//    cyclic_comparator_dynamic0_low1_high2(const node_ptr& header, const node_ptr& x0, node_ptr* path0, node_ptr* path1, node_ptr* path2)
//      : header_(header), x0_(x0), path0_(path0), path1_(path1), path2_(path2) {}
//
//    explicit cyclic_comparator_dynamic0_low1_high2(const node_ptr& header, const node_ptr& x0)
//      : cyclic_comparator_dynamic0_low1_high2(header, x0, default_path_buffer::path0, default_path_buffer::path1, default_path_buffer::path2) { }
//    
//    bool operator()(const node_ptr& x1, const node_ptr& x2) {                  
//      auto iter0 = algo::get_path(x0_, path0_);
//      auto iter1 = algo::get_path(x1, path1_);
//      auto iter2 = algo::get_path(x2, path2_);
//      
//      
//      auto x1Side = x0_ == x1 || algo::less_than(iter0, iter1, path0_);
//      auto x2Side = x0_ == x2 ? false : algo::less_than(iter0, iter2, path0_);
//      return x1Side == x2Side ? algo::less_than(iter1, iter2, path1_) : x1Side;                              
//    }
//  };

//  template<class Algo>
//  class cyclic_comparator_static0_low1_high2
//  {    
//    typedef Algo algo;    
//    typedef typename algo::node_traits node_traits;
//    typedef typename node_traits::node_ptr node_ptr;
//    typedef typename node_traits::node node;
//
//    typedef boost::intrusive::default_path_buffer<node_traits> default_path_buffer;
//
//    node_ptr* path1_;
//    node_ptr* path2_;
//
//  public:    
//    cyclic_comparator_static0_low1_high2(const node_ptr& header, const node_ptr& x0, node_ptr* path1, node_ptr* path2)
//      : path1_(path1), path2_(path2) {                    
//      cyclic_operations<algo>::principal_cut(header, x0);
//    }
//
//    explicit cyclic_comparator_static0_low1_high2(const node_ptr& header, const node_ptr& x0)
//      : cyclic_comparator_static0_low1_high2(header, x0, default_path_buffer::path0, default_path_buffer::path1) { }
//    
//    bool operator()(const node_ptr& x1, const node_ptr& x2) {                  
//      return algo::less_than(algo::get_path(x1, path1_), algo::get_path(x2, path2_), path1_);
//    }
//  };


}}