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

namespace HW::dynamic_connectivity {
  struct etnte_node;
  struct directed_edge;
  struct vertex;
}


namespace HW::dynamic_connectivity
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
    friend struct avl_traits<et_node>;


    /**
     *  has to be the first field for the smart_pool
     *
     *  size_t has to be 8 bytes = 64 bits on a x64 system
     *
     *  [61 bits, element pointer]
     *  [1 bit] - '0' for directed edge, '1' for vertex
     *  [2 bits for avl balance]
     */
    std::enable_if_t<sizeof(size_t) == 8, size_t> ptr_type_balance; 

    // ReSharper disable CppInconsistentNaming
    et_node* parent_;
    et_node* left_;
    et_node* right_;
    // ReSharper restore CppInconsistentNaming

    static inline constexpr auto ptr_bits = ~((static_cast<size_t>(1) << 3) - 1);
    static inline constexpr auto type_bits = static_cast<size_t>(1) << 2;
    static inline constexpr auto balance_bits = type_bits - 1;
  };
  

  


  struct et_traits : avl_traits<et_node>
  {
  private:
    template<class T> 
    static T* get_ptr(const_node_ptr n) {
      return reinterpret_cast<T*>(n->ptr_type_balance & et_node::ptr_bits);  // NOLINT(performance-no-int-to-ptr)
    }

  public:

    // -----------

    static etnte_node* get_non_tree_edge_header(const_node_ptr n) { return get_ptr<etnte_node>(n); }

    static void set_non_tree_edge_header(node_ptr n, etnte_node* etnte_header) {
      n->ptr_type_balance = reinterpret_cast<size_t>(etnte_header) | et_node::balance_bits; // TODO?
    }

    // -----------

    static directed_edge* get_directed_edge(const_node_ptr n) { return get_ptr<directed_edge>(n); }
    static vertex* get_vertex(const_node_ptr n) { return get_ptr<vertex>(n); }
    
    static void set_directed_edge(node_ptr n, directed_edge* directed_edge) {                  
      n->ptr_type_balance = reinterpret_cast<size_t>(directed_edge) | (n->ptr_type_balance & et_node::balance_bits);
    }
    
    static void set_vertex(node_ptr n, vertex* vertex) {
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

    static balance get_balance(const_node_ptr n) {
      return static_cast<balance>(n->ptr_type_balance & et_node::balance_bits);
    }

    static void set_balance(node_ptr n, balance b) {
      n->ptr_type_balance = (n->ptr_type_balance & ~et_node::balance_bits) | static_cast<size_t>(b);
    }   
  };

  using et_node_ptr = et_traits::node_ptr;
  using et_algo = boost::intrusive::avltree_algorithms_ext<et_traits>;









      
  // typedef bi::avl_extended_tree_algorithms<et_traits> euler_tour_algorithms;
  // typedef cyclic_operations<euler_tour_algorithms> euler_tour_cyclic_operations;
  //
  // typedef euler_tour_node_ptr et_node_ptr;
  
  // typedef euler_tour_cyclic_operations et_cyclic_op;
}
