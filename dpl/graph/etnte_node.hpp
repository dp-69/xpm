#pragma once

#include "et_node.hpp"
#include "_OLD_general.hpp"
#include "aug_avltree_algorithms_ext.hpp"
#include "cyclic_operations.hpp"


namespace HW::dynamic_connectivity {
  struct directed_edge;
  using directed_edge_ptr = directed_edge*;
  struct vertex;
}

namespace HW::dynamic_connectivity
{
  struct etnte_traits : dpl::graph::aug_avl_traits<dpl::graph::aug_avl_node>
  {
    static directed_edge_ptr get_directed_edge(const_node_ptr n) {
      return dpl::graph::mask::get_ptr<directed_edge>(n->tag);
    }

    static void set_directed_edge(node_ptr n, const directed_edge_ptr de) {
      dpl::graph::mask::set_ptr(n->tag, de);
    }


    // ordering of non-tree edges
    static et_traits::node_ptr get_vertex_entry(const directed_edge_ptr& de) {
      // sorted by a pointing-in vertex, the pointing-out works as well        
//      return de->opposite_->v1_->et_entry_;   
      return nullptr; // TODO:         return de->v1->et_entry_;
    }
    
    static et_traits::node_ptr get_vertex_entry(const node_ptr& n) {
      return get_vertex_entry(get_directed_edge(n));      
    }
  };


  using etnte_node_ptr = etnte_traits::node_ptr;
  using etnte_algo = boost::intrusive::aug_avltree_algorithms_ext<etnte_traits>;



  
  // inline pair<size_t, size_t> num_vertices_and_edges(et_traits::const_node_ptr header) {    
  //   auto etSize = et_algo::size(header); //et_nt::get_size(et_nt::get_parent(header));    
  //   auto vertexCount = (2 + etSize)/3;    
  //   return make_pair(vertexCount, vertexCount - 1 + etnte_traits::get_size(
  //     etnte_traits::get_parent(et_context_traits::get_non_tree_edge_header(header)))/2);      
  // }
  







  
  struct et_relative_less_than_comparator // key < x
  {
    using et_node_ptr = et_traits::node_ptr;

    et_node_ptr x0_least;
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x1_key;    
    et_traits::const_node_ptr* x1_key_it;

    bool x1_side;

    using default_path = boost::intrusive::default_path_buffer<et_traits>;

    /**
     * x0_least - the least entry, representing a pseudo principal cut
     * x1_key   - is a entry of interest
     * x2       - running entry
     */
    et_relative_less_than_comparator(et_node_ptr least, et_node_ptr key) {
      x0_least = least;
      x0_least_it = et_algo::get_path(x0_least, default_path::path0);

      x1_key = key;
      x1_key_it = et_algo::get_path(x1_key, default_path::path1);      

      x1_side = x0_least == x1_key || et_algo::less_than(x0_least_it, x1_key_it, default_path::path0);
    }

    bool key_less_than_node(const etnte_node_ptr& x2_etnte) const {
      auto x2 = etnte_traits::get_vertex_entry(x2_etnte);

      if (x1_key == x2)
        return false;

      auto x2_it = et_algo::get_path(x2, default_path::path2);

      auto x2_side = x0_least == x2 || et_algo::less_than(x0_least_it, x2_it, default_path::path0);
      return x1_side == x2_side ? et_algo::less_than(x1_key_it, x2_it, default_path::path1) : x1_side;                                    
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






  struct et_relative_more_than_comparator // x < key
  {
    using et_node_ptr = et_traits::node_ptr;

    et_node_ptr x0_least;    
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x2_key;    
    et_traits::const_node_ptr* x2_key_it;

    bool x2_side;

    using default_path = boost::intrusive::default_path_buffer<et_traits>;

    /**
     * x0_least - the least entry, representing a pseudo principal cut
     * x1       - running entry
     * x2_key   - is a entry of interest
     */
    explicit et_relative_more_than_comparator(const et_node_ptr& least, const et_node_ptr& key) {
      x0_least = least;
      x0_least_it = et_algo::get_path(x0_least, default_path::path0);

      x2_key = key;
      x2_key_it = et_algo::get_path(x2_key, default_path::path2);

      x2_side = x0_least == x2_key || et_algo::less_than(x0_least_it, x2_key_it, default_path::path0);
    }

    bool key_more_than_node(const etnte_node_ptr& x1_etnte) const {
      auto x1 = etnte_traits::get_vertex_entry(x1_etnte);

      if (x1 == x2_key)
        return false;

      auto x1_it = et_algo::get_path(x1, default_path::path1);

      auto x1_side = x0_least == x1 || et_algo::less_than(x0_least_it, x1_it, default_path::path0);
      return x1_side == x2_side ? et_algo::less_than(x1_it, x2_key_it, default_path::path1) : x1_side;                                    
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
  



  struct etnte_et_operations
  {
    using et_node_ptr = et_traits::node_ptr;

    using less_than = et_relative_less_than_comparator;
    using more_than = et_relative_more_than_comparator;

    static et_node_ptr get_least_et_entry(const etnte_node_ptr header) {
      return etnte_traits::get_vertex_entry(etnte_traits::get_left(header));
    }

    static etnte_node_ptr lower_bound(etnte_node_ptr header, et_node_ptr vertex_et_entry) {
      return etnte_algo::lower_bound(header, more_than(get_least_et_entry(header), vertex_et_entry));
    }

    static void insert(etnte_node_ptr header, etnte_node_ptr inserting_node) {      
      if (etnte_traits::get_parent(header))
        etnte_algo::insert_equal_upper_bound(header, inserting_node,
          less_than(get_least_et_entry(header), etnte_traits::get_vertex_entry(inserting_node)));

//      Equivalent
//      algo::insert_equal_lower_bound(header, inserting_node,
//        more_than_comparator(get_least_et_entry(header), node_traits::get_vertex_entry(inserting_node)));
      else
        etnte_algo::push_back(header, inserting_node);
    }


    static void split(etnte_node_ptr header_a, etnte_node_ptr header_b, et_node_ptr et_ab, et_node_ptr et_ba) {                
      if (!etnte_traits::get_parent(header_a))
        return;
      
      etnte_node_ptr etnte_least = etnte_traits::get_left(header_a);            
      et_node_ptr et_least = etnte_traits::get_vertex_entry(etnte_least);      

      etnte_node_ptr etnte_entry_ab = etnte_algo::upper_bound(header_a, less_than(et_least, et_ab));
      etnte_node_ptr etnte_entry_ba = etnte_algo::upper_bound(header_a, less_than(et_least, et_ba));     

//      Equivalent
//      auto etnteEntryAB = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ab));
//      auto etnteEntryBA = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ba));     

      
      if (etnte_entry_ab == header_a)
        etnte_entry_ab = etnte_least;

      if (etnte_entry_ba == header_a)
        etnte_entry_ba = etnte_least;

      if (etnte_entry_ab != etnte_entry_ba) {
        if (etnte_algo::less_than(etnte_entry_ab, etnte_entry_ba)) {          
          cyclic_operations<etnte_algo>::split(header_a, header_b, etnte_entry_ab, etnte_entry_ba);
          
          etnte_algo::push_front(header_b, etnte_entry_ab);
          etnte_algo::push_front(header_a, etnte_entry_ba);
        }
        else {
          cyclic_operations<etnte_algo>::split(header_a, header_b, etnte_entry_ba, etnte_entry_ab);
          
          etnte_algo::push_front(header_b, etnte_entry_ba);
          etnte_algo::push_front(header_a, etnte_entry_ab);

          etnte_algo::swap_tree(header_a, header_b);
        }
      }
      else {
        if (cyclic_operations<et_algo>::less_than_low_low(et_least, et_ab, et_ba)) {
          // B is empty                             
        }
        else 
          etnte_algo::swap_tree(header_a, header_b); // A is empty
      }
    }
  };
}
