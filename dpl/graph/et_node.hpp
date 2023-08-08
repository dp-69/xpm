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

#include "avl_extended_augmented_tree_algorithms.hpp"

#include "cyclic_operations.hpp"
#include "dynamic_connectivity_graph.hpp"
#include "smart_pool.hpp"

namespace dpl::graph
{
  using namespace HW;
  using namespace HW::dynamic_connectivity;

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


  struct euler_tour_non_tree_edge_node;
      
  struct et_traits //: default_avl_lrpb_node_traits<et_node>
  {    
  private:
    template<class T> 
    static T* get_ptr(const et_node* n) {
      return reinterpret_cast<T*>(n->ptr_type_balance & et_node::ptr_bits);  // NOLINT(performance-no-int-to-ptr)
    }

  public:
    using node = et_node;
    using node_ptr = node*;
    using const_node_ptr = const node*;

    using balance = avl_balance;

    static constexpr balance negative() { return avl_balance::negative_t; }
    static constexpr balance zero() { return avl_balance::zero_t; }
    static constexpr balance positive() { return avl_balance::positive_t; }

    typedef euler_tour_non_tree_edge_node* etnte_node_ptr;

    // -----------

    static etnte_node_ptr get_non_tree_edge_header(const et_node* n) {
      return get_ptr<euler_tour_non_tree_edge_node>(n);
    }

    static void set_non_tree_edge_header(et_node* n, etnte_node_ptr etnte_header) {
      n->ptr_type_balance = reinterpret_cast<size_t>(etnte_header) | et_node::balance_bits; // TODO?
    }

    // -----------

    static bool is_loop_edge(const et_node* n) { // that is 'vertex'
      return n->ptr_type_balance & et_node::type_bits;
    }
    
    static directed_edge* get_directed_edge(const et_node* n) {
      return get_ptr<directed_edge>(n);
    }

    static void set_directed_edge(et_node* n, directed_edge* directed_edge) {                  
      n->ptr_type_balance = reinterpret_cast<size_t>(directed_edge) | (n->ptr_type_balance & et_node::balance_bits);
    }

    static vertex* get_vertex(const et_node* n) {
      return get_ptr<vertex>(n);
    } 
    
    static void set_vertex(et_node* n, vertex* vertex) {
      n->ptr_type_balance = reinterpret_cast<size_t>(vertex) | et_node::type_bits | (n->ptr_type_balance & et_node::balance_bits);      
    }                     

    // -----------

    static void init_header(et_node* n) {
      n->ptr_type_balance = (n->ptr_type_balance & ~et_node::balance_bits) | et_node::balance_bits; // TODO?
    }

    // Has to be very efficient O(1) predicate.
    static bool is_header(const et_node* n) {
      return (n->ptr_type_balance & et_node::balance_bits) == et_node::balance_bits;
    }

    // -----------

    static node* get_left(const node* n) { return n->left; }
    static node* get_right(const node* n) { return n->right; }
    static node* get_parent(const node* n) { return n->parent; }
    static balance get_balance(const et_node* n) {
      return static_cast<balance>(n->ptr_type_balance & et_node::balance_bits);
    }

    static void set_left(node* n, node* l) { n->left = l; }
    static void set_right(node* n, node* r) { n->right = r; }
    static void set_parent(node* n, node* p) { n->parent = p; }
    static void set_balance(et_node* n, balance b) {
      n->ptr_type_balance = (n->ptr_type_balance & ~et_node::balance_bits) | static_cast<size_t>(b);
    }   
  };

      
  typedef bi::avl_extended_tree_algorithms<et_traits> euler_tour_algorithms;
  typedef cyclic_operations<euler_tour_algorithms> euler_tour_cyclic_operations;
  
  typedef euler_tour_node_ptr et_node_ptr;
  typedef euler_tour_algorithms et_algo;
  typedef euler_tour_cyclic_operations et_cyclic_op;
  




  //  typedef bi::avl_extended_tree_algorithms<euler_tour_node_traits> et_algo_not_augmented;
  



  //  template<>
  //  struct smart_pool_traits<et_node>
  //  {
  //    static et_node_ptr get_next(et_node_ptr x) {            
  //      return static_cast<et_node_ptr>(x->graph_element_entry_type_avl_balance);
  //    }
  //
  //    static void set_next(et_node_ptr x, et_node_ptr y) {
  //      x->graph_element_entry_type_avl_balance = y;
  //    }
  //  };

  //  class euler_tour_iterator : public iterator<forward_iterator_tag, euler_tour_node>
  //  {    
  //    typedef euler_tour_node_traits::node_ptr node_ptr;
  //
  //    node_ptr _node;      
  //
  //  public:
  //    explicit euler_tour_iterator(const node_ptr& node)
  //      : _node(node) {}
  //
  //    euler_tour_iterator& operator++() {
  //      _node = euler_tour_cyclic_operations::next_node(_node);        
  //      return *this;
  //    }
  //
  //    euler_tour_iterator& operator++(int) {
  //      auto temp(*this);
  //      operator++();
  //      return temp;
  //    }
  //
  //    bool operator==(const euler_tour_iterator& rhs) const {
  //      return _node == rhs._node;
  //    }
  //
  //    bool operator!=(const euler_tour_iterator& rhs) const {
  //      return _node != rhs._node;
  //    }
  //
  //    node_ptr ptr() const {
  //      return _node;
  //    }
  //
  //    euler_tour_node& operator*() const {
  //      return *_node;
  //    }
  //
  //    node_ptr operator->() const {
  //      return _node;
  //    }
  //  };
}
