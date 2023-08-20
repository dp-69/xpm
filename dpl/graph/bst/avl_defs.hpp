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

namespace dpl::graph
{
  struct avl_node
  {
    avl_node() = default;
    avl_node(const avl_node& other) = delete;
    avl_node(avl_node&& other) noexcept = delete;
    avl_node& operator=(const avl_node& other) = delete;
    avl_node& operator=(avl_node&& other) noexcept = delete;

    /**
     *  has to be the first field for the smart_pool
     *
     *  must be 8 bytes (64 bits)
     *  61 bits - context data, e.g. pointer
     *  1 bit   - boolean flag
     *  2 bits  - avl balance
     */
    std::enable_if_t<sizeof(std::size_t) == 8, std::size_t> tag; 

    avl_node* parent;
    avl_node* left;
    avl_node* right;
  };


  struct aug_avl_node
  {
    aug_avl_node() = default;
    aug_avl_node(const aug_avl_node& other) = delete;
    aug_avl_node(aug_avl_node&& other) noexcept = delete;
    aug_avl_node& operator=(const aug_avl_node& other) = delete;
    aug_avl_node& operator=(aug_avl_node&& other) noexcept = delete;

    /**
     *  has to be the first field for the smart_pool
     *
     *  must be 8 bytes (64 bits)
     *  61 bits - context data, e.g. pointer
     *  1 bit   - boolean flag
     *  2 bits  - avl balance
     */
    std::enable_if_t<sizeof(std::size_t) == 8, std::size_t> tag;

    aug_avl_node* parent;
    aug_avl_node* left;
    aug_avl_node* right;

    std::size_t size;
  };


  struct mask_bit
  {
    static inline constexpr std::size_t ptr = ~static_cast<std::size_t>(7);
    static inline constexpr std::size_t bit = 4;

    template <typename Ptr> 
    static Ptr get_ptr(std::size_t tag) {
      return reinterpret_cast<Ptr>(tag & ptr);  // NOLINT(performance-no-int-to-ptr)
    }

    template <typename T> 
    static void set_ptr(std::size_t& tag, const T* pointer) {
      tag = reinterpret_cast<size_t>(pointer);
    }

    template <typename T> 
    static void set_ptr_bit(std::size_t& tag, const T* p) {
      tag = reinterpret_cast<size_t>(p) | bit;
    }

    static bool get_bit(std::size_t tag) {
      return tag & bit;
    }
  };


  struct mask_bit_balance : mask_bit
  {
    static inline constexpr std::size_t balance = 3;

    static std::size_t get_balance(std::size_t tag) {
      return tag & balance;
    }

    static void set_balance(std::size_t& tag, std::size_t b) {
      tag = (tag & ~balance) | b;
    }

    template <typename T> 
    static void set_ptr(std::size_t& tag, const T* p) {
      tag = reinterpret_cast<size_t>(p) | (tag & balance);
    }

    template <typename T> 
    static void set_ptr_bit(std::size_t& tag, const T* p) {
      tag = reinterpret_cast<size_t>(p) | bit | (tag & balance);
    }

    template <typename T> 
    static void set_ptr_balance(std::size_t& tag, const T* p, std::size_t b) {
      tag = reinterpret_cast<size_t>(p) | b;
    }
  };


  template<typename Node>
  class avl_traits 
  {
    using mask = mask_bit_balance;

    enum class avl_balance
    {
      negative_t = 0,
      zero_t = 1,
      positive_t = 2,
      fourth_state_t = 3
    };

  public:
    using node = Node;
    using node_ptr = node*;
    using const_node_ptr = const node*;

    using balance = avl_balance;

    static node_ptr get_left(const_node_ptr n) { return n->left; }
    static node_ptr get_right(const_node_ptr n) { return n->right; }
    static node_ptr get_parent(const_node_ptr n) { return n->parent; }

    static void set_left(node_ptr n, node_ptr l) { n->left = l; }
    static void set_right(node_ptr n, node_ptr r) { n->right = r; }
    static void set_parent(node_ptr n, node_ptr p) { n->parent = p; }

    static balance negative() { return balance::negative_t; }
    static balance zero() { return balance::zero_t; }
    static balance positive() { return balance::positive_t; }

    static balance get_balance(const_node_ptr n) {
      return static_cast<balance>(mask::get_balance(n->tag));
    }

    static void set_balance(node_ptr n, balance b) {
      return mask::set_balance(n->tag, static_cast<size_t>(b));
    }

    /**
     * \brief avl balance state of a header equals to (4 in decimal) or (11 as bitmap), that is not a valid avl balance state
     */
    static void init_header(node_ptr n) {
      mask::set_balance(n->tag, mask::balance);
    }

    /**
     * \brief checks balance bits, when equal to 4 (not a valid avl balance state) the node is a header
     */
    static bool is_header(const_node_ptr n) {
      return mask::get_balance(n->tag) == mask::balance;
    }
  };


  template<typename Node>
  struct aug_avl_traits : avl_traits<Node>
  {
    using subsize = std::size_t;

    static subsize get_size(typename avl_traits<Node>::const_node_ptr n) {
      return n->size;
    }

    static void set_size(typename avl_traits<Node>::node_ptr n, subsize s) {
      n->size = s;
    }
  };
}