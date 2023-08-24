﻿#pragma once

#include "et_node.hpp"
#include "_OLD_general.hpp"
#include "avl_extended_augmented_tree_algorithms.hpp"
#include "cyclic_operations.hpp"


namespace HW::dynamic_connectivity {
  struct directed_edge;
  using directed_edge_ptr = directed_edge*;
  struct vertex;
}

namespace HW { namespace dynamic_connectivity
{
  struct etnte_node : non_copyable_movable
  {    
    private:
      friend struct etnte_traits;      
      friend struct avl_traits<etnte_node>;

      // ReSharper disable CppInconsistentNaming


      /**
       *  has to be the first field for the smart_pool
       *
       *  size_t has to be 8 bytes = 64 bits on a x64 system
       *
       *  [61 bits, element pointer]
       *  [1 bit, EMPTY]
       *  [2 bits for avl balance]
       */
      
      std::size_t de_and_balance_;

      etnte_node* parent_;
      etnte_node* left_;
      etnte_node* right_;

      std::size_t size_;
      // ReSharper restore CppInconsistentNaming
  };
  


  // typedef euler_tour_non_tree_edge_node* euler_tour_non_tree_edge_node_ptr;
  
  
  struct etnte_traits : avl_traits<etnte_node>
  {
    using node = etnte_node;
    using node_ptr = etnte_node*;
    using const_node_ptr = const etnte_node*;            
    
    typedef tagged_pointer_as_size_t<balance, 2> compression;

    static void init_header(const node_ptr& n) {
      n->de_and_balance_ = 0;            
    }

    static bool is_header(const const_node_ptr& n) {
      return !n->de_and_balance_;            
    }

    
    static balance get_balance(const const_node_ptr& n) {      
      return compression::get_bits(n->de_and_balance_);
    }

    static void set_balance(const node_ptr& n, balance b) {
      compression::set_bits(n->de_and_balance_, b);
    }

    static directed_edge_ptr get_directed_edge(const node_ptr& n) {      
      return compression::get_pointer<directed_edge_ptr>(n->de_and_balance_);
    }

    static void set_directed_edge(const node_ptr& n, const directed_edge_ptr& de) {     
      compression::set_pointer(n->de_and_balance_, de);      
    }


    // ordering of non-tree edges
    static et_node_ptr get_vertex_entry(const directed_edge_ptr& de) {
      // sorted by a pointing-in vertex, the pointing-out works as well        
//      return de->opposite_->v1_->et_entry_;   
      return nullptr; // TODO:         return de->v1->et_entry_;
    }
    
    static et_node_ptr get_vertex_entry(const node_ptr& n) {
      return get_vertex_entry(get_directed_edge(n));      
    }


    using subsize = size_t;

    static subsize get_size(const_node_ptr n) { return n->size_; /*return n ? n->size_ : 0;*/ }
    static void set_size(node_ptr n, subsize s) { n->size_ = s; }
  };


  using etnte_node_ptr = etnte_traits::node_ptr;
  using etnte_algo = boost::intrusive::avl_extended_augmented_tree_algorithms<etnte_traits>;



  
  inline pair<size_t, size_t> num_vertices_and_edges(et_traits::const_node_ptr header) {    
    auto etSize = et_algo::size(header); //et_nt::get_size(et_nt::get_parent(header));    
    auto vertexCount = (2 + etSize)/3;    
    return make_pair(vertexCount, vertexCount - 1 + etnte_traits::get_size(
      etnte_traits::get_parent(et_traits::get_non_tree_edge_header(header)))/2);      
  }
  






  
  template <class NodeTraits>
  struct euler_tour_pseudo_cut_less_than_comparator // key < x
  {
    et_node_ptr x0_;
    et_traits::const_node_ptr* x0_it_;

    et_node_ptr x1_;    
    et_traits::const_node_ptr* x1_it_;

    bool x1Side_;

    using default_path = boost::intrusive::default_path_buffer<et_traits>;

    // x0_least - the least entry, representing a pseudo principal cut
    // x1_key - is a entry of interest
    // x2/node - running entry

    explicit euler_tour_pseudo_cut_less_than_comparator(const et_node_ptr& x0_least, const et_node_ptr& x1_key) {
      x0_ = x0_least;
      x0_it_ = et_algo::get_path(x0_, default_path::path0);

      x1_ = x1_key;
      x1_it_ = et_algo::get_path(x1_, default_path::path1);      

      x1Side_ = x0_ == x1_ || et_algo::less_than(x0_it_, x1_it_, default_path::path0);
    }

    bool key_less_than_node(const etnte_node_ptr& x2_etnte) const {
      auto x2_ = NodeTraits::get_vertex_entry(x2_etnte);

      if (x1_ == x2_)
        return false;

      auto x2_it = et_algo::get_path(x2_, default_path::path2);

      auto x2Side = x0_ == x2_ || et_algo::less_than(x0_it_, x2_it, default_path::path0);
      return x1Side_ == x2Side ? et_algo::less_than(x1_it_, x2_it, default_path::path1) : x1Side_;                                    
    }

    bool operator()(const etnte_node_ptr& x2_etnte) const {
      return key_less_than_node(x2_etnte);
    }

    bool operator()(const et_node_ptr& /*comparing key*/, const etnte_node_ptr& x2_etnte) const {
      return key_less_than_node(x2_etnte);
    }

    bool operator()(const etnte_node_ptr& /*comparing node*/, const etnte_node_ptr& x2_etnte) const {
      return key_less_than_node(x2_etnte);    
    }
  }; 






  template <class NodeTraits>
  struct euler_tour_pseudo_cut_more_than_comparator // x < key
  {
    et_node_ptr x0_;    
    et_traits::const_node_ptr* x0_it_;

    et_node_ptr x2_;    
    et_traits::const_node_ptr* x2_it_;

    bool x2Side_;

    using default_path = boost::intrusive::default_path_buffer<et_traits>;

    // x0_least - the least entry, representing a pseudo principal cut
    // x1/node - running entry
    // x2_key - is a entry of interest

    explicit euler_tour_pseudo_cut_more_than_comparator(const et_node_ptr& x0_least, const et_node_ptr& x2_key) {
      x0_ = x0_least;
      x0_it_ = et_algo::get_path(x0_, default_path::path0);

      x2_ = x2_key;
      x2_it_ = et_algo::get_path(x2_, default_path::path2);

      x2Side_ = x0_ == x2_ || et_algo::less_than(x0_it_, x2_it_, default_path::path0);
    }

    bool key_more_than_node(const etnte_node_ptr& x1_etnte) const {
      auto x1_ = NodeTraits::get_vertex_entry(x1_etnte);

      if (x1_ == x2_)
        return false;

      auto x1_it = et_algo::get_path(x1_, default_path::path1);

      auto x1Side = x0_ == x1_ || et_algo::less_than(x0_it_, x1_it, default_path::path0);
      return x1Side == x2Side_ ? et_algo::less_than(x1_it, x2_it_, default_path::path1) : x1Side;                                    
    }

    bool operator()(const etnte_node_ptr& x1_etnte) const {
      return key_more_than_node(x1_etnte);
    }

    bool operator()(const etnte_node_ptr& x1_etnte, const et_node_ptr& /*comparing key*/) const {
      return key_more_than_node(x1_etnte);
    }

    bool operator()(const etnte_node_ptr& x1_etnte, const etnte_node_ptr& /*comparing node*/) const {
      return key_more_than_node(x1_etnte);    
    }
  };
  



  template<class Algo>
  struct euler_tour_non_tree_edge_operations
  {
    typedef Algo algo;
    using node_traits = typename algo::node_traits;
    using node = typename node_traits::node;
    using node_ptr = typename node_traits::node_ptr;
    using const_node_ptr = typename node_traits::const_node_ptr;
    
    typedef euler_tour_pseudo_cut_less_than_comparator<node_traits> less_than_comparator;
    typedef euler_tour_pseudo_cut_more_than_comparator<node_traits> more_than_comparator;
    
    typedef cyclic_operations<algo> cyclic_op;

    static et_node_ptr get_least_et_entry(const node_ptr& header) {
      return node_traits::get_vertex_entry(node_traits::get_left(header));
    }

    static node_ptr lower_bound(const node_ptr& header, const et_node_ptr& vertex_et_entry) {
      return algo::lower_bound(header, more_than_comparator(get_least_et_entry(header), vertex_et_entry));
    }

    static void insert(const node_ptr& header, const node_ptr& inserting_node) {      
      if (node_traits::get_parent(header))
        algo::insert_equal_upper_bound(header, inserting_node,
          less_than_comparator(get_least_et_entry(header), node_traits::get_vertex_entry(inserting_node)));

//      Equivalent
//      algo::insert_equal_lower_bound(header, inserting_node,
//        more_than_comparator(get_least_et_entry(header), node_traits::get_vertex_entry(inserting_node)));
      else
        algo::push_back(header, inserting_node);
    }


    static void split(const node_ptr& headerA, const node_ptr& headerB, et_node_ptr et_ab, et_node_ptr et_ba) {                
      if (!etnte_traits::get_parent(headerA))
        return;

      auto etnteLeast = node_traits::get_left(headerA);            
      auto etLeast = node_traits::get_vertex_entry(etnteLeast);      

      auto etnteEntryAB = algo::upper_bound(headerA, less_than_comparator(etLeast, et_ab));
      auto etnteEntryBA = algo::upper_bound(headerA, less_than_comparator(etLeast, et_ba));     

//      Equivalent
//      auto etnteEntryAB = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ab));
//      auto etnteEntryBA = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ba));     

      
      if (etnteEntryAB == headerA)
        etnteEntryAB = etnteLeast;

      if (etnteEntryBA == headerA)
        etnteEntryBA = etnteLeast;

      if (etnteEntryAB != etnteEntryBA) {
        if (algo::less_than(etnteEntryAB, etnteEntryBA)) {          
          cyclic_op::split(headerA, headerB, etnteEntryAB, etnteEntryBA);
          
          algo::push_front(headerB, etnteEntryAB);
          algo::push_front(headerA, etnteEntryBA);
        }
        else {
          cyclic_op::split(headerA, headerB, etnteEntryBA, etnteEntryAB);
          
          algo::push_front(headerB, etnteEntryBA);
          algo::push_front(headerA, etnteEntryAB);

          algo::swap_tree(headerA, headerB);
        }
      }
      else {
        if (cyclic_operations<et_algo>::less_than_low_low(etLeast, et_ab, et_ba)) {
          // B is empty                             
        }
        else 
          algo::swap_tree(headerA, headerB); // A is empty
      }
    }

    
  };

  
  
  
  
  
  typedef euler_tour_non_tree_edge_operations<etnte_algo> etnte_op;
  typedef cyclic_operations<etnte_algo> etnte_cyclic_op;
}}
