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

#include "core.hpp"
#include "avltree_algorithms_ext.hpp"

namespace HW::dynamic_connectivity
{
  struct et_node
  {
    et_node() = default;
    et_node(const et_node& other) = delete;
    et_node(et_node&& other) noexcept = delete;
    et_node& operator=(const et_node& other) = delete;
    et_node& operator=(et_node&& other) noexcept = delete;

    /**
     *  has to be the first field for the smart_pool
     *
     *  must be 8 bytes (64 bits)
     *  61 bits - context data, e.g. pointer
     *  1 bit   - boolean flag
     *  2 bits  - avl balance
     */
    std::enable_if_t<sizeof(std::size_t) == 8, std::size_t> tag; 

    et_node* parent;
    et_node* left;
    et_node* right;
  };

  using et_traits = avl_traits<et_node>;
  using et_node_ptr = et_traits::node_ptr;
  using et_algo = boost::intrusive::avltree_algorithms_ext<et_traits>;
}









namespace HW::dynamic_connectivity
{
  struct etnte_node;
  struct directed_edge;
  struct vertex;

  class et_context_traits
  {
    using node_ptr = et_node*;
    using const_node_ptr = const et_node*;

  public:
    static etnte_node* get_non_tree_edge_header(const_node_ptr n) {
      return mask::get_ptr<etnte_node>(n->tag);
    }

    static void set_non_tree_edge_header(node_ptr n, etnte_node* etnte_header) {
      mask::set_ptr_balance(n->tag, etnte_header, mask::balance);
    }

    static directed_edge* get_directed_edge(const_node_ptr n) {
      return mask::get_ptr<directed_edge>(n->tag);
    }

    static vertex* get_vertex(const_node_ptr n) {
      return mask::get_ptr<vertex>(n->tag);
    }
    
    static void set_directed_edge(node_ptr n, const directed_edge* de) {
      mask::set_ptr(n->tag, de);
    }
    
    static void set_vertex(node_ptr n, const vertex* v) {
      mask::set_ptr_bit(n->tag, v);
    }                     

    static bool is_loop_edge(const_node_ptr n) { // that is 'vertex'
      return mask::get_bit(n->tag);
    }
  };
}
