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

namespace HW::dynamic_connectivity {
  struct euler_tour_non_tree_edge_node;
  struct directed_edge;
  struct vertex;
}


namespace dpl::graph
{
  struct et_node
  {
    et_node() = default;
    et_node(const et_node& other) = delete;
    et_node(et_node&& other) noexcept = delete;
    et_node& operator=(const et_node& other) = delete;
    et_node& operator=(et_node&& other) noexcept = delete;

  private:
    friend struct et_traits;

    // has to be the first field for the smart_pool
    // size_t has to be 8 bytes = 64 bits, x64 system
    // [61 bits, element pointer] + [1 bit, entry type] + [2 bits for avl balance]
    std::enable_if_t<sizeof(size_t) == 8, size_t> ptr_type_balance; 

    et_node* parent;
    et_node* left;
    et_node* right;

    static inline constexpr auto ptr_bits = ~((static_cast<size_t>(1) << 3) - 1);
    static inline constexpr auto type_bits = static_cast<size_t>(1) << 2;
    static inline constexpr auto balance_bits = type_bits - 1;
  };

  


  struct et_traits
  {
    using node = et_node;
    using node_ptr = node*;
    using const_node_ptr = const node*;

    using balance = avl_balance;

    static constexpr balance negative() { return avl_balance::negative_t; }
    static constexpr balance zero() { return avl_balance::zero_t; }
    static constexpr balance positive() { return avl_balance::positive_t; }

  private:
    template<class T> 
    static T* get_ptr(const_node_ptr n) {
      return reinterpret_cast<T*>(n->ptr_type_balance & et_node::ptr_bits);  // NOLINT(performance-no-int-to-ptr)
    }

    using de = HW::dynamic_connectivity::directed_edge;
    using ve = HW::dynamic_connectivity::vertex;
    using pt = HW::dynamic_connectivity::euler_tour_non_tree_edge_node;

  public:

    // -----------

    static pt* get_non_tree_edge_header(const_node_ptr n) { return get_ptr<pt>(n); }

    static void set_non_tree_edge_header(node_ptr n, pt* etnte_header) {
      n->ptr_type_balance = reinterpret_cast<size_t>(etnte_header) | et_node::balance_bits; // TODO?
    }

    // -----------

    static de* get_directed_edge(const_node_ptr n) { return get_ptr<de>(n); }
    static ve* get_vertex(const_node_ptr n) { return get_ptr<ve>(n); }
    
    static void set_directed_edge(node_ptr n, de* directed_edge) {                  
      n->ptr_type_balance = reinterpret_cast<size_t>(directed_edge) | (n->ptr_type_balance & et_node::balance_bits);
    }
    
    static void set_vertex(node_ptr n, ve* vertex) {
      n->ptr_type_balance = reinterpret_cast<size_t>(vertex) | et_node::type_bits | (n->ptr_type_balance & et_node::balance_bits);      
    }                     


    static bool is_loop_edge(const_node_ptr n) { // that is 'vertex'
      return n->ptr_type_balance & et_node::type_bits;
    }

    // -----------

    static void init_header(node_ptr n) {
      n->ptr_type_balance = (n->ptr_type_balance & ~et_node::balance_bits) | et_node::balance_bits; // TODO?
    }

    // Has to be very efficient O(1) predicate.
    static bool is_header(const_node_ptr n) {
      return (n->ptr_type_balance & et_node::balance_bits) == et_node::balance_bits;
    }

    // -----------

    static node_ptr get_left(const_node_ptr n) { return n->left; }
    static node_ptr get_right(const_node_ptr n) { return n->right; }
    static node_ptr get_parent(const_node_ptr n) { return n->parent; }
    static balance get_balance(const_node_ptr n) {
      return static_cast<balance>(n->ptr_type_balance & et_node::balance_bits);
    }

    static void set_left(node_ptr n, node_ptr l) { n->left = l; }
    static void set_right(node_ptr n, node_ptr r) { n->right = r; }
    static void set_parent(node_ptr n, node_ptr p) { n->parent = p; }
    static void set_balance(node_ptr n, balance b) {
      n->ptr_type_balance = (n->ptr_type_balance & ~et_node::balance_bits) | static_cast<size_t>(b);
    }   
  };

      
  // typedef bi::avl_extended_tree_algorithms<et_traits> euler_tour_algorithms;
  // typedef cyclic_operations<euler_tour_algorithms> euler_tour_cyclic_operations;
  //
  // typedef euler_tour_node_ptr et_node_ptr;
  // typedef euler_tour_algorithms et_algo;
  // typedef euler_tour_cyclic_operations et_cyclic_op;
}
