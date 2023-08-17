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

#include "bst/avl_defs.hpp"
#include "bst/aug_avltree_algorithms_ext.hpp"
#include "bst/cyclic.hpp"

namespace HW::dynamic_connectivity
{
  using et_traits = dpl::graph::avl_traits<dpl::graph::avl_node>;
  using et_algo = dpl::graph::avltree_algorithms_ext<et_traits>;


  using etnte_traits = dpl::graph::aug_avl_traits<dpl::graph::aug_avl_node>;
  using etnte_algo = dpl::graph::aug_avltree_algorithms_ext<etnte_traits>;


  template <typename Context>
  struct et_relative_less_than_comparator // key < x
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    et_node_ptr x0_least;
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x1_key;    
    et_traits::const_node_ptr* x1_key_it;

    bool x1_side;

    using default_path = dpl::graph::default_path_buffer<et_traits>;

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
      auto x2 = Context::get_vertex_entry(x2_etnte);

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


  template <typename Context>
  struct et_relative_more_than_comparator // x < key
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    et_node_ptr x0_least;    
    et_traits::const_node_ptr* x0_least_it;

    et_node_ptr x2_key;    
    et_traits::const_node_ptr* x2_key_it;

    bool x2_side;

    using default_path = dpl::graph::default_path_buffer<et_traits>;

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
      auto x1 = Context::get_vertex_entry(x1_etnte);

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


  template <typename Context>
  struct etnte_context_operations
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    using less_than = et_relative_less_than_comparator<Context>;
    using more_than = et_relative_more_than_comparator<Context>;

    static et_node_ptr get_least_et_entry(const etnte_node_ptr header) {
      return Context::get_vertex_entry(etnte_traits::get_left(header));
    }

    static etnte_node_ptr lower_bound(etnte_node_ptr header, et_node_ptr vertex_et_entry) {
      return etnte_algo::lower_bound(header, more_than(get_least_et_entry(header), vertex_et_entry));
    }

    static void insert(etnte_node_ptr header, etnte_node_ptr inserting_node) {      
      if (etnte_traits::get_parent(header))
        etnte_algo::insert_equal_upper_bound(header, inserting_node,
          less_than(get_least_et_entry(header), Context::get_vertex_entry(inserting_node)));

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
      et_node_ptr et_least = Context::get_vertex_entry(etnte_least);      

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
          dpl::graph::cyclic<etnte_algo>::split(header_a, header_b, etnte_entry_ab, etnte_entry_ba);
          
          etnte_algo::push_front(header_b, etnte_entry_ab);
          etnte_algo::push_front(header_a, etnte_entry_ba);
        }
        else {
          dpl::graph::cyclic<etnte_algo>::split(header_a, header_b, etnte_entry_ba, etnte_entry_ab);
          
          etnte_algo::push_front(header_b, etnte_entry_ba);
          etnte_algo::push_front(header_a, etnte_entry_ab);

          etnte_algo::swap_tree(header_a, header_b);
        }
      }
      else {
        if (dpl::graph::cyclic<et_algo>::less_than_low_low(et_least, et_ab, et_ba)) {
          // B is empty                             
        }
        else 
          etnte_algo::swap_tree(header_a, header_b); // A is empty
      }
    }
  };
}


 
  // inline pair<size_t, size_t> num_vertices_and_edges(et_traits::const_node_ptr header) {    
  //   auto etSize = et_algo::size(header); //et_nt::get_size(et_nt::get_parent(header));    
  //   auto vertexCount = (2 + etSize)/3;    
  //   return make_pair(vertexCount, vertexCount - 1 + etnte_traits::get_size(
  //     etnte_traits::get_parent(et_context_traits::get_non_tree_edge_header(header)))/2);      
  // }