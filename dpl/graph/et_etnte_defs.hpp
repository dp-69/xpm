/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * DPL is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * DPL is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "bst/aug_avltree_algorithms_ext.hpp"
#include "bst/avl_defs.hpp"
#include "bst/cyclic.hpp"

#include <boost/graph/depth_first_search.hpp>
#include <boost/pool/object_pool.hpp>

namespace dpl::graph::helper
{
  template <class T,
    typename = std::enable_if_t<sizeof(T) >= sizeof(std::size_t)>>
  class smart_pool
  { 
    static T* get_next(T* x) {       
      return *reinterpret_cast<T**>(x);
    }

    static void set_next(T* x, T* next) {
      *reinterpret_cast<T**>(x) = next;
    }

    boost::object_pool<T> pool_;
    T* released_stack_top_ = nullptr;

  public:
    T* acquire() {
      if (released_stack_top_) {
        auto result = released_stack_top_;        
        released_stack_top_ = get_next(released_stack_top_);       
        return result;
      }      

      return pool_.malloc();
    }

    void release(T* x) {
      set_next(x, released_stack_top_);
      released_stack_top_ = x;
    }
  };
}

namespace dpl::graph
{
  using et_traits = avl_traits<avl_node>;
  using et_algo = avltree_algorithms_ext<et_traits>;


  using etnte_traits = aug_avl_traits<aug_avl_node>;
  using etnte_algo = aug_avltree_algorithms_ext<etnte_traits>;



  template <typename Props>
  void print(etnte_traits::node_ptr hdr, const Props& c) {
    for (etnte_traits::node_ptr etnte : range<etnte_traits>(hdr)) {
      auto de = Props::get_directed_edge(etnte);
      std::cout << fmt::format(" ({}, {})", c.get_idx(Props::get_opposite(de)->v1), c.get_idx(de->v1));
    }
  }

  template <typename Props>
  void print(et_traits::node_ptr hdr, const Props& c) {
    for (et_traits::node_ptr et : range<et_traits>(hdr)) {
      if (Props::is_loop_edge(et)) {
        auto v = Props::get_vertex(et);
        std::cout << fmt::format(" [{}]", c.get_idx(v));  
      }
      else {
        auto de = Props::get_directed_edge(et);
        std::cout << fmt::format(" ({}, {})", c.get_idx(Props::get_opposite(de)->v1), c.get_idx(de->v1));  
      }
    }
  }


  template <typename Props>
  struct et_relative_less_than_comparator // key < x
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    et_node_ptr x0_least;
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x1_key;    
    et_traits::const_node_ptr* x1_key_it;

    bool x1_side;

    using default_path = default_path_buffer<et_traits>;

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
      auto x2 = Props::get_ordering_vertex_entry(x2_etnte);

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


  template <typename Props>
  struct et_relative_more_than_comparator // x < key
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    et_node_ptr x0_least;    
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x2_key;    
    et_traits::const_node_ptr* x2_key_it;

    bool x2_side;

    using default_path = default_path_buffer<et_traits>;

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
      auto x1 = Props::get_ordering_vertex_entry(x1_etnte);

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


  template <typename Props>
  struct etnte_context_operations
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    using less_than = et_relative_less_than_comparator<Props>;
    using more_than = et_relative_more_than_comparator<Props>;

    static etnte_node_ptr lower_bound(etnte_node_ptr hdr, et_node_ptr vertex_et_entry) {
      return etnte_algo::lower_bound(hdr, more_than{
        Props::get_ordering_vertex_entry(etnte_traits::get_left(hdr)),
        vertex_et_entry
      });
    }

    static void insert(etnte_node_ptr hdr, etnte_node_ptr inserting_node) {      
      if (etnte_traits::get_parent(hdr))
        etnte_algo::insert_equal_upper_bound(hdr, inserting_node,
          less_than{
            Props::get_ordering_vertex_entry(etnte_traits::get_left(hdr)),
            Props::get_ordering_vertex_entry(inserting_node)
          });

//      Equivalent
//      algo::insert_equal_lower_bound(header, inserting_node,
//        more_than_comparator(get_least_et_entry(header), node_traits::get_vertex_entry(inserting_node)));
      else
        etnte_algo::push_back(hdr, inserting_node);
    }


    static void split(etnte_node_ptr header_a, etnte_node_ptr header_b, et_node_ptr et_ab, et_node_ptr et_ba) {                
      if (!etnte_traits::get_parent(header_a))
        return;
      
      etnte_node_ptr etnte_least = etnte_traits::get_left(header_a);            
      et_node_ptr et_least = Props::get_ordering_vertex_entry(etnte_least);      

      etnte_node_ptr etnte_entry_ab = etnte_algo::upper_bound(header_a, less_than{et_least, et_ab});
      etnte_node_ptr etnte_entry_ba = etnte_algo::upper_bound(header_a, less_than{et_least, et_ba});     

//      Equivalent
//      auto etnteEntryAB = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ab));
//      auto etnteEntryBA = algo::lower_bound(headerA, more_than_comparator(etLeast, et_ba));     

      
      if (etnte_entry_ab == header_a)
        etnte_entry_ab = etnte_least;

      if (etnte_entry_ba == header_a)
        etnte_entry_ba = etnte_least;

      if (etnte_entry_ab != etnte_entry_ba) {
        if (etnte_algo::less_than(etnte_entry_ab, etnte_entry_ba)) {
          cyclic<etnte_algo>::split(header_a, header_b, etnte_entry_ab, etnte_entry_ba);
          
          etnte_algo::push_front(header_b, etnte_entry_ab);
          etnte_algo::push_front(header_a, etnte_entry_ba);
        }
        else {
          cyclic<etnte_algo>::split(header_a, header_b, etnte_entry_ba, etnte_entry_ab);
          
          etnte_algo::push_front(header_b, etnte_entry_ba);
          etnte_algo::push_front(header_a, etnte_entry_ab);

          etnte_algo::swap_tree(header_a, header_b);
        }
      }
      else {
        if (cyclic<et_algo>::less_than_low_low(et_least, et_ab, et_ba)) {
          // B is empty                             
        }
        else 
          etnte_algo::swap_tree(header_a, header_b); // A is empty
      }
    }
  };


  template <typename Graph, typename Props>
  class euler_tour_visitor : public boost::default_dfs_visitor
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    using v_desc = typename boost::graph_traits<Graph>::vertex_descriptor;
    using e_desc = typename boost::graph_traits<Graph>::edge_descriptor;

    et_node_ptr et_hdr_;
    etnte_node_ptr etnte_hdr_;

    helper::smart_pool<et_traits::node>* et_pool_;
    helper::smart_pool<etnte_traits::node>* etnte_pool_;

    et_node_ptr* tree_edge_stack_empty_;
    et_node_ptr* tree_edge_stack_top_;
    
  public:
    euler_tour_visitor(
      helper::smart_pool<et_traits::node>* et_pool
    , helper::smart_pool<etnte_traits::node>* etnte_pool
    , et_node_ptr* tree_edge_stack
    )
      : et_pool_(et_pool)
      , etnte_pool_(etnte_pool)
      , tree_edge_stack_empty_(tree_edge_stack)
      , tree_edge_stack_top_(tree_edge_stack)
    {}


    void start_vertex(v_desc, const Graph&) {            
      et_hdr_ = et_pool_->acquire();
      etnte_hdr_ = etnte_pool_->acquire();
      et_algo::init_header(et_hdr_);      
      etnte_algo::init_header(etnte_hdr_);      
      Props::set_etnte_header(et_hdr_, etnte_hdr_);             
    }

    void discover_vertex(v_desc v, const Graph&) {     
      et_node_ptr entry = et_pool_->acquire();        
      Props::set_vertex(entry, v);
      Props::set_entry(v, entry);
      et_algo::push_back(et_hdr_, entry);            
    }

    void finish_vertex(v_desc v, const Graph&) {
      if (tree_edge_stack_top_ != tree_edge_stack_empty_) {
        et_node_ptr entry = et_pool_->acquire();
        et_node_ptr top = *tree_edge_stack_top_--;
        e_desc top_de_opposite = Props::get_opposite(Props::get_directed_edge(top));
        Props::set_directed_edge(entry, top_de_opposite);        
        Props::set_tree_edge_entry(top_de_opposite, entry);
        et_algo::push_back(et_hdr_, entry);        
      }
    }    

    void tree_edge(e_desc e, const Graph&) {
      et_node_ptr entry = et_pool_->acquire();
      Props::set_directed_edge(entry, e);
      Props::set_tree_edge_entry(e, entry);
      et_algo::push_back(et_hdr_, entry);      
      *++tree_edge_stack_top_ = entry;
    }

    // non-tree edge    
    void forward_or_cross_edge(e_desc de, const Graph&) {
      etnte_node_ptr entry = etnte_pool_->acquire();
      etnte_node_ptr entry_opp = etnte_pool_->acquire();
      
      e_desc de_opp = Props::get_opposite(de);

      Props::set_non_tree_edge_entry(de, entry); 
      Props::set_non_tree_edge_entry(de_opp, entry_opp);
      
      Props::set_directed_edge(entry, de);
      Props::set_directed_edge(entry_opp, de_opp);
            
      etnte_context_operations<Props>::insert(etnte_hdr_, entry);
      etnte_context_operations<Props>::insert(etnte_hdr_, entry_opp);                    
    }
  }; 
}


 
  // inline pair<size_t, size_t> num_vertices_and_edges(et_traits::const_node_ptr header) {    
  //   auto etSize = et_algo::size(header); //et_nt::get_size(et_nt::get_parent(header));    
  //   auto vertexCount = (2 + etSize)/3;    
  //   return make_pair(vertexCount, vertexCount - 1 + etnte_traits::get_size(
  //     etnte_traits::get_parent(et_context_traits::get_non_tree_edge_header(header)))/2);      
  // }