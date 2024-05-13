/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "dc_traits.hpp"
#include "bst/cyclic.hpp"

#include <boost/graph/two_bit_color_map.hpp>
#include <boost/property_map/function_property_map.hpp>

namespace dpl::graph
{      
  // struct component_visitor_empty
  // {
  //   void vertex(vertex*) {}
  //   void directed_edge(directed_edge*) {}
  // };
    
  // struct component_visitor_count
  // {
  //   size_t _vertices = 0;
  //   size_t _directedEdges = 0;
  //
  //   void vertex(vertex*) { ++_vertices; }
  //   void directed_edge(directed_edge*) { ++_directedEdges; }
  //
  //   size_t vertices() { return _vertices; }
  //   size_t edges() { return _directedEdges/2; }
  // };


  // struct dc_stat
  // {
  //   size_t nteCountOther;
  //   size_t nteCountActive;
  //   size_t nteReconnectingFirst;
  // };

  // class nte_iterator
  // {
  //   dc_graph* graph_;
  //
  //   using et_iter_t = internal::bst_inorder_iterator<et_traits>;
  //   et_iter_t et_iter, et_end;
  //
  //   out_edge_iterator edge_iter, edge_end;
  //
  // public:
  //   using edge_t = dc_graph::edge_descriptor;
  //   
  //   nte_iterator() = default;
  //   nte_iterator(et_traits::node_ptr hdr, dc_graph* g) {
  //     graph_ = g;
  //
  //     et_iter = et_iter_t{et_traits::get_left(hdr)};
  //     et_end = et_iter_t{hdr};
  //
  //     while (true) {
  //       if (end())
  //         return;
  //
  //       if (dc_properties::is_loop_edge(*et_iter)) {
  //         auto v_idx = dc_properties::get_vertex(*et_iter);
  //
  //         for (edge_iter = graph_->out_edges_begin(v_idx), edge_end = graph_->out_edges_end(v_idx); edge_iter != edge_end; ++edge_iter)
  //           if (!dc_properties::is_null_entry(*edge_iter, *graph_) && !dc_properties::is_tree_edge(*edge_iter, *graph_))
  //             return;
  //       }
  //
  //       ++et_iter;
  //     }
  //   }
  //
  //   edge_t operator*() const { return *edge_iter; }
  //   edge_t operator->() const { return *edge_iter; }
  //
  //   nte_iterator& operator++() {
  //     ++edge_iter;
  //
  //     for (; edge_iter != edge_end; ++edge_iter)
  //       if (!dc_properties::is_null_entry(*edge_iter, *graph_) && !dc_properties::is_tree_edge(*edge_iter, *graph_))
  //         return *this;
  //
  //     ++et_iter;
  //
  //     while (true) {
  //       if (end())
  //         return *this;
  //
  //       if (dc_properties::is_loop_edge(*et_iter)) {
  //         auto v_idx = dc_properties::get_vertex(*et_iter);
  //
  //         for (edge_iter = graph_->out_edges_begin(v_idx), edge_end = graph_->out_edges_end(v_idx); edge_iter != edge_end; ++edge_iter)
  //           if (!dc_properties::is_null_entry(*edge_iter, *graph_) && !dc_properties::is_tree_edge(*edge_iter, *graph_))
  //             return *this;
  //       }
  //
  //       ++et_iter;
  //     }
  //   }
  //
  //   bool end() const {
  //     return et_iter == et_end;
  //   }
  // };

  class dc_context
  {
    using et_ptr = et_traits::node_ptr;
    using et_nt = et_traits;

    using vertex_t = dc_graph::vertex_t;
    using edge_t = dc_graph::edge_t;

    dc_graph* g_;

    static inline constexpr auto null_edge_ = std::numeric_limits<edge_t>::max();

  public:

    helper::smart_pool<et_traits::node> et_pool_;

    void init_with_dfs(dc_graph& g) {
      g_ = &g;

      auto vertex_count = g.vertex_count();

      auto color_uptr = std::make_unique<boost::two_bit_color_type[]>(vertex_count);

      auto color_map = boost::make_function_property_map<dc_graph::vertex_t>(
        [&](dc_graph::vertex_t v) -> boost::two_bit_color_type& { return color_uptr[*v]; });

      auto tree_edge_stack = std::vector<et_ptr>(vertex_count + 1);     
      euler_tour_visitor visitor{&et_pool_, tree_edge_stack.data()};

      for (dc_graph::vertex_t i{0}; i < g.vertex_count(); ++i) {
        if (color_map[i] != boost::two_bit_color_type::two_bit_black) {
          boost::depth_first_visit(g, i, visitor, color_map); // tree edges

          for (et_ptr et : range<et_nt>(et_algo::get_header(get_entry(i, *g_))))
            if (is_loop_edge(et, *g_))
              for (edge_t de : g.edges(get_vertex(et, *g_)))
                if (is_null_entry(de, *g_)) // non-tree edges
                  set_non_tree_edge(de, *g_);
        }
      }
    }
    
    void adjacent_edges_remove(vertex_t v, const dc_graph& g) {
      for (edge_t de : g.edges(v))
        if (!is_null_entry(de, g) && !is_tree_edge(de, g))
          non_tree_edge_remove(de);

      for (edge_t de : g.edges(v))
        if (!is_null_entry(de, g))
          tree_edge_split_and_reconnect(de);

      set_entry(v, nullptr, g);
      // TODO: release vertex?
    }

    void non_tree_edge_remove(edge_t ab) {
      set_null_entry(ab, *g_);
      set_null_entry(opposite(ab, *g_), *g_);
    }    


    edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
      if (
        et_algo::node_height(et_nt::get_parent(hdr_a)) <
        et_algo::node_height(et_nt::get_parent(hdr_b)))
        for (et_ptr et_entry : range<et_nt>(hdr_a)) {
          if (is_loop_edge(et_entry, *g_)) {
            for (edge_t ab : g_->edges(get_vertex(et_entry, *g_)))
              if (!is_null_entry(ab, *g_) && !is_tree_edge(ab, *g_))
                if (et_algo::get_header(get_entry(target(ab, *g_), *g_)) == hdr_b)
                  return ab;
          }
        }
      else
        for (et_ptr et_entry : range<et_nt>(hdr_b)) {
          if (is_loop_edge(et_entry, *g_)) {
            for (edge_t ba : g_->edges(get_vertex(et_entry, *g_)))
              if (!is_null_entry(ba, *g_) && !is_tree_edge(ba, *g_))
                if (et_algo::get_header(get_entry(target(ba, *g_), *g_)) == hdr_a)
                  return opposite(ba, *g_);
          }
        }
    
      return null_edge_;
    }



    // edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
    //   auto [iter_a, end_a] = range<et_nt>(hdr_a);
    //   auto [iter_b, end_b] = range<et_nt>(hdr_b);
    //
    //   while (true) {
    //     if (iter_a == end_a || iter_b == end_b)
    //       return nullptr;
    //
    //     if (dc_properties::is_loop_edge(*iter_a)) {
    //       for (edge_t ab : graph_->edges(dc_properties::get_vertex(*iter_a)))
    //         if (!dc_properties::is_null_entry(ab) && !dc_properties::is_tree_edge(ab))
    //           if (et_algo::get_header(Props::get_entry(ab->v1)) == hdr_b)
    //             return ab;
    //     }
    //
    //     if (dc_properties::is_loop_edge(*iter_b)) {
    //       for (edge_t ba : graph_->edges(dc_properties::get_vertex(*iter_b)))
    //         if (!dc_properties::is_null_entry(ba) && !dc_properties::is_tree_edge(ba))
    //           if (et_algo::get_header(Props::get_entry(ba->v1)) == hdr_a)
    //             return Props::get_opposite(ba);
    //     }
    //
    //     ++iter_a;
    //     ++iter_b;
    //   }
    // }

    


    // edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
    //   nte_iterator iter_ab{hdr_a, graph_};
    //   nte_iterator iter_ba{hdr_b, graph_};
    //
    //   while (true) {
    //     if (iter_ab.end() || iter_ba.end())
    //       return nullptr;
    //
    //     if (et_algo::get_header(Props::get_entry(iter_ab->v1)) == hdr_b)
    //       return *iter_ab;
    //
    //     if (et_algo::get_header(Props::get_entry(iter_ba->v1)) == hdr_a)
    //       return Props::get_opposite(*iter_ba);
    //
    //     ++iter_ab;
    //     ++iter_ba;
    //   }
    // }

    /**
     * \return true if reconnected
     */
    bool tree_edge_split_and_reconnect(edge_t ab) { // TODO: check release of entries
      auto ba = opposite(ab, *g_);

      et_ptr et_ab = get_tree_edge_entry(ab, *g_);
      et_ptr et_ba = get_tree_edge_entry(ba, *g_);
      
      et_ptr et_hdr_a = et_algo::get_header(et_ab);
      et_ptr et_hdr_b = et_pool_.acquire();
      et_algo::init_header(et_hdr_b);

      if (et_algo::less_than(et_ab, et_ba))
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ab, et_ba);
      else {
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ba, et_ab);
        et_algo::swap_tree(et_hdr_a, et_hdr_b);        
      }      
                                                                
      set_null_entry(ab, *g_);
      set_null_entry(ba, *g_);





      // if (auto quick_idx = saved_replacement_indices[quick_replacement_idx++];
      //   quick_idx != std::numeric_limits<std::size_t>::max())

      // replace_ab_entry = 
      //   Props::get_non_tree_edge_entry(
      //     graph_->get_directed_edge(quick_idx));


      // for (etnte_ptr etnte : range<etnte_traits>(etnte_hdr_a))
      //   if (etnte_algo::get_header(
      //     Props::get_non_tree_edge_entry(
      //       Props::get_opposite(
      //         Props::get_directed_edge(etnte)))) == etnte_hdr_b) {
      //     
      //     replace_ab_entry = etnte;
      //     break;
      //   }



      // saved_replacement_indices.push_back(
      //   replace_ab_entry
      //     ? graph_->get_idx(
      //       swap
      //         ? Props::get_directed_edge(replace_ab_entry)->opposite
      //         : Props::get_directed_edge(replace_ab_entry)
      //       ) 
      //     : std::numeric_limits<std::size_t>::max()
      // );

      


      if (ab = find_replacement(et_hdr_a, et_hdr_b); ab != null_edge_) {                                
        ba = opposite(ab, *g_);        

        et_ptr replace_a_entry = get_ordering_vertex_entry(ab, *g_);
        et_ptr replace_b_entry = get_ordering_vertex_entry(ba, *g_);                                        

        set_directed_edge(et_ab, ab, *g_);
        set_directed_edge(et_ba, ba, *g_);

        set_tree_edge_entry(ab, et_ab, *g_);
        set_tree_edge_entry(ba, et_ba, *g_);

        // et
        cyclic<et_algo>::cut(et_hdr_a, replace_a_entry);
        cyclic<et_algo>::cut(et_hdr_b, replace_b_entry);                
                  
        et_algo::join_trees(et_hdr_a, et_ab, et_hdr_b);        
        et_algo::push_back(et_hdr_a, et_ba);                       

        et_pool_.release(et_hdr_b);

        return true;
      }

      et_pool_.release(et_ab);
      et_pool_.release(et_ba);

      return false;
    }


        


    void release_component(et_ptr hdr) { // TODO: check consistency
      et_ptr et_node = et_nt::get_left(hdr);      
      while (et_node != hdr) {
        if (is_loop_edge(et_node, *g_))
          set_entry(get_vertex(et_node, *g_), nullptr, *g_);
        else
          set_null_entry(get_directed_edge(et_node, *g_), *g_);

        et_node = et_algo::next_node(et_node);
      }

      et_node = et_nt::get_left(hdr);
      while (et_node != hdr) {
        auto prev = et_node;
        et_node = et_algo::next_node(et_node);
        et_pool_.release(prev);
      }

      et_pool_.release(hdr);
    }
  };  
}
