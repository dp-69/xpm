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
#include "dc_graph.hpp"

#include <boost/graph/depth_first_search.hpp>
#include <boost/pool/object_pool.hpp>

namespace dpl::graph
{
  namespace helper
  {
    template <class T, typename = std::enable_if_t<sizeof(T) >= sizeof(std::size_t)>>
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

  using et_traits = avl_traits<avl_node>;
  using et_algo = avltree_algorithms_ext<et_traits>;

  inline et_traits::node_ptr get_tree_edge_entry(dc_graph::edge_t x, const dc_graph& g) {
    return mask_bit::get_ptr<et_traits::node_ptr>(g.directed_edge_entry(x));
  }

  inline void set_tree_edge_entry(dc_graph::edge_t x, et_traits::node_ptr y, const dc_graph& g) {
    mask_bit::set_ptr(g.directed_edge_entry(x), y);
  }

  inline void set_non_tree_edge(dc_graph::edge_t x, const dc_graph& g) { // TODO
    return mask_bit::set_bit(g.directed_edge_entry(x));
  }

  inline dc_graph::vertex_t get_vertex(et_traits::const_node_ptr n, const dc_graph&) {
    return dc_graph::vertex_t{mask_bit_balance::get_value<dc_graph::vertex_t::value_type>(n->tag)};
  }

  inline void set_vertex(et_traits::node_ptr n, dc_graph::vertex_t v, const dc_graph&) {
    mask_bit_balance::set_value_bit(n->tag, *v);
  }

  inline dc_graph::edge_t get_directed_edge(et_traits::const_node_ptr n, const dc_graph&) {
    return dc_graph::edge_t{mask_bit_balance::get_value<dc_graph::edge_t::value_type>(n->tag)};
  }

  inline void set_directed_edge(et_traits::node_ptr n, dc_graph::edge_t de, const dc_graph&) {
    mask_bit_balance::set_value(n->tag, *de);
  }

  /**
   * \brief checks if n is a vertex-type et entry
   */
  inline bool is_loop_edge(et_traits::const_node_ptr n, const dc_graph&) {
    return mask_bit_balance::get_bit(n->tag);
  }

  inline void set_entry(dc_graph::vertex_t v, et_traits::node_ptr et, const dc_graph& g) {
    g.vertex_entry(v) = et;
  }

  inline et_traits::node_ptr get_entry(dc_graph::vertex_t v, const dc_graph& g) {
    return static_cast<et_traits::node_ptr>(g.vertex_entry(v));
  }

  /**
   * \brief
   *   ordering of non-tree edges
   *   sorted by a pointing-in vertex, the pointing-out works as well
   */
  inline et_traits::node_ptr get_ordering_vertex_entry(dc_graph::edge_t de, const dc_graph& g) {
    return get_entry(target(opposite(de, g), g), g);
  }

  inline bool is_tree_edge(dc_graph::edge_t x, const dc_graph& g) {
    return !mask_bit::get_bit(g.directed_edge_entry(x));
  }

  inline bool is_null_entry(dc_graph::edge_t x, const dc_graph& g) {
    return g.directed_edge_entry(x) == 0;
  }

  inline void set_null_entry(dc_graph::edge_t x, const dc_graph& g) {
    g.directed_edge_entry(x) = 0;
  }






  class euler_tour_visitor : public boost::default_dfs_visitor
  {
    using et_ptr = et_traits::node_ptr;

    using vertex_t = dc_graph::vertex_t;
    using edge_t = dc_graph::edge_t;

    et_ptr et_hdr_;

    helper::smart_pool<et_traits::node>* et_pool_;

    et_ptr* tree_edge_stack_empty_;
    et_ptr* tree_edge_stack_top_;
    
  public:
    euler_tour_visitor(
      helper::smart_pool<et_traits::node>* et_pool,
      et_ptr* tree_edge_stack
    ) : et_pool_(et_pool),
        tree_edge_stack_empty_(tree_edge_stack),
        tree_edge_stack_top_(tree_edge_stack)
    {}


    void start_vertex(vertex_t, const dc_graph&) {            
      et_hdr_ = et_pool_->acquire();
      et_algo::init_header(et_hdr_);      
    }

    void discover_vertex(vertex_t v, const dc_graph& g) {     
      et_ptr et = et_pool_->acquire();        
      set_vertex(et, v, g);
      set_entry(v, et, g);
      et_algo::push_back(et_hdr_, et);            
    }

    void finish_vertex(vertex_t, const dc_graph& g) {
      if (tree_edge_stack_top_ != tree_edge_stack_empty_) {
        et_ptr top = *tree_edge_stack_top_--;
        et_ptr et = et_pool_->acquire();
        edge_t top_de_opposite = opposite(get_directed_edge(top, g), g);
        set_directed_edge(et, top_de_opposite, g);        
        set_tree_edge_entry(top_de_opposite, et, g);
        et_algo::push_back(et_hdr_, et);
      }
    }    

    void tree_edge(edge_t e, const dc_graph& g) {
      et_ptr et = et_pool_->acquire();
      set_directed_edge(et, e, g);
      set_tree_edge_entry(e, et, g);
      et_algo::push_back(et_hdr_, et);      
      *++tree_edge_stack_top_ = et;
    }
  }; 
}