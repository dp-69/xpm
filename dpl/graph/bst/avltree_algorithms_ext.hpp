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

#include <boost/intrusive/avltree_algorithms.hpp>

namespace dpl::graph
{
  template<typename NodeTraits>
  struct default_path_buffer // TODO: maybe move to the function iself, check performance
  {
    using node_ptr = typename NodeTraits::node_ptr;
    using const_node_ptr = typename NodeTraits::const_node_ptr;

    static inline constexpr auto ET_MAX_DEPTH = 256;

    static inline const_node_ptr path0[ET_MAX_DEPTH];
    static inline const_node_ptr path1[ET_MAX_DEPTH];
    static inline const_node_ptr path2[ET_MAX_DEPTH];
    static inline const_node_ptr path3[ET_MAX_DEPTH];
  };


  template <typename NodeTraits>
  class avltree_algorithms_ext : public boost::intrusive::avltree_algorithms<NodeTraits>
  {
    using base = boost::intrusive::avltree_algorithms<NodeTraits>;
    using nt = NodeTraits;

  public:
    using node_traits = NodeTraits;

    using node = typename nt::node;
    using node_ptr = typename nt::node_ptr;
    using const_node_ptr = typename nt::const_node_ptr;
    using balance = typename nt::balance;
    
    static node_ptr get_header(const_node_ptr node) {
      while (!is_header(node))
        node = nt::get_parent(node);
      return const_cast<node_ptr>(node);
    }

    static void init(node_ptr node) {
      base::init(node);
      nt::set_balance(node, nt::zero());
    }

    static void init_header(node_ptr header) {
      boost::intrusive::bstree_algorithms<nt>::init_header(header);        
      nt::init_header(header);
    }

    static bool is_header(const_node_ptr p) {
      return nt::is_header(p);
    }

    static bool empty(const_node_ptr header) {
      return !nt::get_parent(header);
    }

    static auto node_height(const_node_ptr x) {
      int height = 0;

      while (nt::get_left(x) || nt::get_right(x)) {
        x = nt::get_balance(x) == nt::negative()
              ? nt::get_left(x)
              : nt::get_right(x);

        ++height;
      }

      return height;
    }

    static const_node_ptr* get_path(const_node_ptr n, const_node_ptr* path) {   
      while (!is_header(n)) {
        *++path = n; // TODO
        n = nt::get_parent(n);
      }       

      return path;
    }      
    
    /**
     * \brief l and r must not be equal
     */
    static bool less_than(const_node_ptr l, const_node_ptr r) {
      using default_path = default_path_buffer<nt>;

      auto l_subpath = get_path(l, default_path::path0);
      auto r_subpath = get_path(r, default_path::path1);
      return less_than(l_subpath, r_subpath, default_path::path0);
    }

    static bool less_than(const_node_ptr* l_subpath, const_node_ptr* r_subpath, const_node_ptr* l_path) {  // TODO                      
      // equal if *(path0 + 1) == *(path1 + 1) // TODO - try to understand
      // should not be equal
      
      while (*l_subpath == *r_subpath) {
        --l_subpath;
        --r_subpath;
      }

      if (l_subpath == l_path)
        return nt::get_right(*(r_subpath + 1)) == *r_subpath;

      return nt::get_left(*(l_subpath + 1)) == *l_subpath;        
    }    

    template<typename UnaryLessThanNodeComparator>
    static node_ptr upper_bound(node_ptr header, UnaryLessThanNodeComparator comp) {
      node_ptr x = nt::get_parent(header);
      node_ptr y = header;

      while (x)
        if (comp(x)) {     // Key < x  
          y = x;
          x = nt::get_left(x);
        }
        else
          x = nt::get_right(x);

      return y;         
    }

    template<typename UnaryMoreThanNodeComparator>
    static node_ptr lower_bound(node_ptr header, UnaryMoreThanNodeComparator comp) {
      node_ptr x = nt::get_parent(header);
      node_ptr y = header;

      while (x)
        if (comp(x))   // x < Key 
          x = nt::get_right(x);
        else {
          y = x;
          x = nt::get_left(x);
        }

      return y;         
    }

    static size_t calculate_subtree_size(const_node_ptr p) {
      return base::subtree_size(p);
    }

    static void join_trees(node_ptr hdr_left, node_ptr x, node_ptr hdr_right) {
      node_ptr root_a = nt::get_parent(hdr_left);
      node_ptr root_b = nt::get_parent(hdr_right);

      if (!root_b) {
        base::push_back(hdr_left, x);
        return;
      }

      if (!root_a) {
        base::swap_tree(hdr_left, hdr_right);            
        base::push_front(hdr_left, x);
        return;
      }

      auto height_a = node_height(root_a);
      auto height_b = node_height(root_b);            
        
      if (std::abs(height_b - height_a) <= 1) {                        
        nt::set_parent(hdr_left, x);
        nt::set_left(hdr_left, nt::get_left(hdr_left));
        nt::set_right(hdr_left, nt::get_right(hdr_right));

        nt::set_parent(x, hdr_left);

        nt::set_left(x, root_a);
        nt::set_right(x, root_b);

        nt::set_parent(root_a, x);
        nt::set_parent(root_b, x);

        nt::set_balance(x,
          height_a < height_b ? nt::positive() :
          height_a > height_b ? nt::negative() :
          nt::zero()
        );

        init_header(hdr_right);
        return;            
      }
        

      if (height_a > height_b) {
        node_ptr v = root_a;
        auto h = height_a;
        while (h > height_b + 1) {
          if (nt::get_balance(v) == nt::negative())
            h -= 2;
          else
            --h;
          v = nt::get_right(v);
        }

        node_ptr v_parent = nt::get_parent(v);

        nt::set_left(x, v);
        nt::set_parent(v, x);

        nt::set_right(x, root_b);        
        nt::set_parent(root_b, x);

        nt::set_right(v_parent, x);
        nt::set_parent(x, v_parent);

        nt::set_balance(x, height_b == h ? nt::zero() : nt::negative());

        nt::set_right(hdr_left, nt::get_right(hdr_right));
  
        init_header(hdr_right);
      }
      else {            
        node_ptr v = root_b;
        auto h = height_b;
        while (h > height_a + 1) {
          if (nt::get_balance(v) == nt::positive())
            h -= 2;
          else
            --h;

          v = nt::get_left(v);
        }

        node_ptr v_parent = nt::get_parent(v);

        nt::set_right(x, v);
        nt::set_parent(v, x);

        nt::set_left(x, root_a);        
        nt::set_parent(root_a, x);
  
        nt::set_left(v_parent, x);
        nt::set_parent(x, v_parent);

        nt::set_balance(x, height_a == h ? nt::zero() : nt::positive());

        nt::set_left(hdr_right, nt::get_left(hdr_left));
  
        init_header(hdr_left);
        base::swap_tree(hdr_left, hdr_right);
      }

      base::rebalance_after_insertion_no_balance_assignment(hdr_left, x);
    }


    static void split_tree(node_ptr hdr_left, node_ptr k, node_ptr hdr_right) {
      if (nt::get_right(hdr_left) == k) {
        base::erase(hdr_left, k);
        return;
      }

      if (nt::get_left(hdr_left) == k) {
        base::erase(hdr_left, k);
        base::swap_tree(hdr_left, hdr_right);
        return;
      }

      auto k_height = node_height(k);

      node_ptr left_tail_node = nt::get_left(k);
      auto left_tail_height = k_height - (nt::get_balance(k) == nt::positive() ? 2 : 1);

      node_ptr right_tail_node = nt::get_right(k);
      auto right_tail_height = k_height - (nt::get_balance(k) == nt::negative() ? 2 : 1);

      node_ptr left_rightmost = nullptr;
      node_ptr right_leftmost = nullptr;

      if (left_tail_node) {
        left_rightmost = left_tail_node;
        while (nt::get_right(left_rightmost))
          left_rightmost = nt::get_right(left_rightmost);
      }

      if (right_tail_node) {
        right_leftmost = right_tail_node;
        while (nt::get_left(right_leftmost))
          right_leftmost = nt::get_left(right_leftmost);
      }


      node_ptr node = k;
      node_ptr parent = nt::get_parent(k);
      node_ptr grand_parent = nt::get_parent(parent);

      auto height = k_height;

      while (parent != hdr_left) {
        if (nt::get_left(parent) == node) {
          if (!right_leftmost)
            right_leftmost = parent;

          auto l_height = right_tail_height;
          auto r_height =
            nt::get_balance(parent) == nt::positive()
              ? ++++height - 1
              : nt::get_balance(parent) == nt::negative()
                  ? ++height - 2
                  : ++height - 1;

          node_ptr x = parent;
          node_ptr l_node = right_tail_node;
          node_ptr r_node = nt::get_right(x);

          auto height_diff = r_height - l_height;

          if (height_diff <= 1) {
            nt::set_left(x, l_node);
            if (l_node)
              nt::set_parent(l_node, x);

            if (height_diff == 1) {
              nt::set_balance(x, nt::positive());
              right_tail_height = r_height + 1;
            }
            else if (height_diff == 0) {
              nt::set_balance(x, nt::zero());
              right_tail_height = r_height + 1;
            }
            else /*if (heightDiff == -1)*/ {
              nt::set_balance(x, nt::negative());
              right_tail_height = r_height + 2;
            }

            right_tail_node = x;
          }
          else {
            if (l_height == -1) {
              nt::set_left(x, nullptr);
              nt::set_right(x, nullptr);
              nt::set_balance(x, nt::zero());
              
              node_ptr l = r_node;

              while (nt::get_left(l))
                l = nt::get_left(l);

              nt::set_left(l, x);
              nt::set_parent(x, l);
            }
            else {
              node_ptr v = r_node;
              auto h = r_height;
              while (h > l_height + 1) {
                if (nt::get_balance(v) == nt::positive())
                  h -= 2;
                else
                  --h;                                    

                v = nt::get_left(v);
              }

              node_ptr v_parent = nt::get_parent(v);

              nt::set_right(x, v);
              nt::set_parent(v, x);

              nt::set_left(x, l_node);
              if (l_node)
                nt::set_parent(l_node, x);

              nt::set_left(v_parent, x);
              nt::set_parent(x, v_parent);

              nt::set_balance(x, l_height == h ? nt::zero() : nt::positive());
            }

            nt::set_parent(hdr_right, r_node);
            nt::set_parent(r_node, hdr_right);

            if (base::rebalance_after_insertion_no_balance_assignment(hdr_right, x))
              right_tail_height = r_height;
            else
              right_tail_height = r_height + 1;

            right_tail_node = nt::get_parent(hdr_right);
          }
        }
        else {
          if (!left_rightmost)
            left_rightmost = parent;

          auto l_height =
            nt::get_balance(parent) == nt::negative()
              ? ++++height - 1
              : nt::get_balance(parent) == nt::positive()
                  ? ++height - 2
                  : ++height - 1;

          auto r_height = left_tail_height;

          node_ptr x = parent;
          node_ptr l_node = nt::get_left(x);
          node_ptr r_node = left_tail_node;

          auto height_diff = l_height - r_height;

          if (height_diff <= 1) {
            nt::set_right(x, r_node);
            if (r_node)
              nt::set_parent(r_node, x);

            if (height_diff == 1) {
              nt::set_balance(x, nt::negative());
              left_tail_height = l_height + 1;
            }
            else if (height_diff == 0) {
              nt::set_balance(x, nt::zero());
              left_tail_height = l_height + 1;
            }
            else /*if (heightDiff == -1)*/ {
              nt::set_balance(x, nt::positive());
              left_tail_height = l_height + 2;
            }

            left_tail_node = x;
          }
          else {
            if (r_height == -1) {
              nt::set_left(x, nullptr);
              nt::set_right(x, nullptr);
              nt::set_balance(x, nt::zero());
              
              node_ptr r = l_node;

              while (nt::get_right(r))
                r = nt::get_right(r);

              nt::set_right(r, x);
              nt::set_parent(x, r);               
            }
            else {
              node_ptr v = l_node;
              auto h = l_height;
              while (h > r_height + 1) {
                if (nt::get_balance(v) == nt::negative())
                  h -= 2;
                else
                  --h;

                v = nt::get_right(v);
              }

              node_ptr v_parent = nt::get_parent(v);

              nt::set_left(x, v);
              nt::set_parent(v, x);

              nt::set_right(x, r_node);
              if (r_node)
                nt::set_parent(r_node, x);

              nt::set_right(v_parent, x);
              nt::set_parent(x, v_parent);

              nt::set_balance(x, r_height == h ? nt::zero() : nt::negative());
            }

            nt::set_parent(hdr_right, l_node);
            nt::set_parent(l_node, hdr_right);

            if (base::rebalance_after_insertion_no_balance_assignment(hdr_right, x))
              left_tail_height = l_height;
            else
              left_tail_height = l_height + 1;

            left_tail_node = nt::get_parent(hdr_right);
          }
        }

        node = parent;
        parent = grand_parent;
        grand_parent = nt::get_parent(grand_parent);
      }

      nt::set_parent(hdr_right, right_tail_node);
      nt::set_left(hdr_right, right_leftmost);
      nt::set_right(hdr_right, nt::get_right(hdr_left));

      nt::set_parent(right_tail_node, hdr_right);


      nt::set_parent(hdr_left, left_tail_node);
      /* nt::set_left(headerA, nt::get_left(headerA)); */                    
      nt::set_right(hdr_left, left_rightmost);

      nt::set_parent(left_tail_node, hdr_left);

      init(k);
    }
  };
}
