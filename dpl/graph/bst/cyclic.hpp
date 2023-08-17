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
  template <typename Algo>
  class cyclic
  {
    using algo = Algo;
    using nt = typename algo::node_traits;

    using node = typename nt::node;
    using node_ptr = typename nt::node_ptr;
    using const_node_ptr = typename nt::const_node_ptr;

  public:    
    static void principal_cut(node_ptr header, node_ptr least) {      
      if (nt::get_left(header) != least)
        if (nt::get_right(header) == least) {
          algo::erase(header, least);
          algo::push_front(header, least);
        }
        else {  // NOLINT(clang-diagnostic-dangling-else)
          node s;
          node_ptr hdr_s = &s;
          algo::init_header(hdr_s);          
          algo::split_tree(header, least, hdr_s);
          algo::push_front(hdr_s, least);

          node_ptr right_s = nt::get_right(hdr_s);
          algo::erase(hdr_s, right_s);
          algo::join_trees(hdr_s, right_s, header);
          algo::swap_tree(header, hdr_s);
        }
    }

    static void principal_cut_least_dropped(node_ptr header, node_ptr least) {
      if (nt::get_right(header) == least || nt::get_left(header) == least)
        algo::erase(header, least);
      else {        
        node s;
        auto hdr_s = &s;
        algo::init_header(hdr_s);        
        algo::split_tree(header, least, hdr_s);

        auto right_b = nt::get_right(hdr_s);
        algo::erase(hdr_s, right_b);
        algo::join_trees(hdr_s, right_b, header);
        algo::swap_tree(header, hdr_s);
      }  
    }

    static node_ptr next_node(node_ptr node) {
      node_ptr next = algo::next_node(node);      
      return nt::is_header(next) ? nt::get_left(next) : next;       
    }


    // Split initial BST into the following format: A + {n0} + B + {n1} + C.
    // Returns:
    //   headerB <- B (from n0 to n1)
    //   headerA <- C + A (from n1 to n0)

    /**
     * \brief Assuming format A + {n0} + B + {n1} + C
     * \param hdr_a returns (n1, n0) = C + A range
     * \param hdr_b returns (n0, n1) = B range
     */
    static void split(node_ptr hdr_a, node_ptr hdr_b, node_ptr n0, node_ptr n1) {            
      algo::split_tree(hdr_a, n0, hdr_b);
            
      node c;
      node_ptr hdr_c = &c;
      algo::init_header(hdr_c);      
      algo::split_tree(hdr_b, n1, hdr_c);

      if (nt::get_parent(hdr_a)) {
        node_ptr left_a = nt::get_left(hdr_a);
        algo::erase(hdr_a, left_a);
        algo::join_trees(hdr_c, left_a, hdr_a);
      }
      
      algo::swap_tree(hdr_a, hdr_c);      
    }
    

    static bool less_than_low_low(const node_ptr& x0, const node_ptr& x1, const node_ptr& x2) {
      using default_path = default_path_buffer<nt>;

      auto iter0 = algo::get_path(x0, default_path::path0);
      auto iter1 = algo::get_path(x1, default_path::path1);
      auto iter2 = algo::get_path(x2, default_path::path2);

      auto x1_side = x0 == x1 || algo::less_than(iter0, iter1, default_path::path0);
      auto x2_side = x0 == x2 || algo::less_than(iter0, iter2, default_path::path0);
      return x1_side == x2_side ? algo::less_than(iter1, iter2, default_path::path1) : x1_side;    
    }
  };
}
