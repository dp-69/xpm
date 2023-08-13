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
#include "avl_defs.hpp"

namespace HW::dynamic_connectivity
{
  using et_traits = dpl::graph::avl_traits<dpl::graph::avl_node>;
  using et_algo = boost::intrusive::avltree_algorithms_ext<et_traits>;
}

namespace HW::dynamic_connectivity
{
  struct directed_edge;
  struct vertex;

  class et_context_traits
  {
    using et_node_ptr = et_traits::node_ptr;
    using et_const_node_ptr = et_traits::const_node_ptr;

    using etnte_node = dpl::graph::aug_avl_node;

    using mask = dpl::graph::mask;

  public:
    static etnte_node* get_non_tree_edge_header(et_const_node_ptr n) {
      return mask::get_ptr<etnte_node>(n->tag);
    }

    static void set_non_tree_edge_header(et_node_ptr n, const etnte_node* etnte_header) {
      mask::set_ptr_balance(n->tag, etnte_header, mask::balance);
    }

    static directed_edge* get_directed_edge(et_const_node_ptr n) {
      return mask::get_ptr<directed_edge>(n->tag);
    }

    static vertex* get_vertex(et_const_node_ptr n) {
      return mask::get_ptr<vertex>(n->tag);
    }
    
    static void set_directed_edge(et_node_ptr n, const directed_edge* de) {
      mask::set_ptr(n->tag, de);
    }
    
    static void set_vertex(et_node_ptr n, const vertex* v) {
      mask::set_ptr_bit(n->tag, v);
    }                     

    static bool is_loop_edge(et_const_node_ptr n) { // that is 'vertex'
      return mask::get_bit(n->tag);
    }
  };
}
