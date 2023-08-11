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

namespace HW::dynamic_connectivity
{
  enum class avl_balance
  {
    negative_t,
    zero_t,
    positive_t,
    fourth_state_t
  };

  template<class Node>
  struct avl_traits 
  {  
    using node = Node;
    using node_ptr = node*;
    using const_node_ptr = const node*;

    typedef avl_balance balance;

    static node_ptr get_left(const_node_ptr n) { return n->left_; }
    static node_ptr get_right(const_node_ptr n) { return n->right_; }
    static node_ptr get_parent(const_node_ptr n) { return n->parent_; }

    static void set_left(node_ptr n, node_ptr l) { n->left_ = l; }
    static void set_right(node_ptr n, node_ptr r) { n->right_ = r; }
    static void set_parent(node_ptr n, node_ptr p) { n->parent_ = p; }

    static balance negative() { return avl_balance::negative_t; }
    static balance zero() { return avl_balance::zero_t; }
    static balance positive() { return avl_balance::positive_t; }
  };
}