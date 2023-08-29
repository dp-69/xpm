#pragma once

#include "dc_properties.hpp"

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


  struct dc_stat
  {
    size_t nteCountOther;
    size_t nteCountActive;
    size_t nteReconnectingFirst;
  };

  template <typename Props>
  class dc_context
  {
    using et_ptr = et_traits::node_ptr;
    using etnte_ptr = etnte_traits::node_ptr;

    using et_nt = et_traits;
    using etnte_nt = etnte_traits;

    using vertex_t = dc_graph::vertex_descriptor;
    using edge_t = dc_graph::edge_descriptor;

    dc_graph* graph_;
    // std::unique_ptr<std::list<edge_t>[]> nte_edges;

  public:
    // std::vector<std::size_t> saved_replacement_indices;
    // std::size_t quick_replacement_idx = 0;

    helper::smart_pool<et_traits::node> et_pool_;
    helper::smart_pool<etnte_traits::node> etnte_pool_;

    /**
     * \brief dfs for et-only, fast direct etnte
     */
    void init_with_dfs(dc_graph& g, dc_properties c) {
      graph_ = &g;
      // saved_replacement_indices.reserve(g.edge_count());


      auto vertex_count = num_vertices(g);                 

      // nte_edges = std::make_unique<std::list<edge_t>[]>(vertex_count);

      using color_t = boost::two_bit_color_type;
      auto color_uptr = std::make_unique<color_t[]>(vertex_count);
      boost::iterator_property_map color_map{
        color_uptr.get(),
        boost::make_function_property_map<vertex_t>(
          [&c](vertex_t v) { return c.get_idx(v); })
      };

      auto tree_edge_stack = std::vector<et_ptr>(vertex_count + 1);     
      euler_tour_visitor<dc_graph, dc_properties> visitor{&et_pool_, &etnte_pool_, tree_edge_stack.data()};

      for (std::size_t i = 0; i < g.vertex_count(); ++i) {
        // nte_edges[i].clear();
        // for (auto de : g.edges(i)) {
        //   nte_edges[i].push_back(de);
        // }

        vertex_t v = g.get_vertex(i);

        if (color_map[v] != color_t::two_bit_black) {
          boost::depth_first_visit(g, v, visitor, color_map); // tree edges

          et_ptr et_hdr = et_algo::get_header(dc_properties::get_entry(v));
          etnte_ptr etnte_hdr = dc_properties::get_etnte_header(et_hdr);

          for (et_ptr et_entry : range<et_nt>(et_hdr))
            if (dc_properties::is_loop_edge(et_entry))
              for (edge_t de : g.edges(dc_properties::get_vertex(et_entry)))
                if (dc_properties::is_null_entry(de)) { // non-tree edges
                  etnte_ptr etnte_entry = etnte_pool_.acquire();
                  dc_properties::set_non_tree_edge_entry(de, etnte_entry); 
                  // dc_properties::set_directed_edge(etnte_entry, de); //TODO

                  // etnte_algo::push_back(etnte_hdr, etnte_entry);
                  // boost::intrusive::avltree_algorithms<etnte_traits>::push_back(etnte_hdr, etnte_entry);


                }

          #ifndef ETNTE_AS_AVL_ONLY
          etnte_algo::populate_sizes(etnte_hdr);
          #endif
        }
      }
    }

    void adjacent_edges_remove(std::integral auto v_idx, const dc_graph& g) {
      for (edge_t de : g.edges(v_idx))
        if (!dc_properties::is_null_entry(de) && !dc_properties::is_tree_edge(de))
          non_tree_edge_remove(de);

      for (edge_t de : g.edges(v_idx))
        if (!dc_properties::is_null_entry(de))
          tree_edge_split_and_reconnect(de);

      // TODO: release vertex?
    }

    void non_tree_edge_remove(edge_t ab) {
      edge_t ba = Props::get_opposite(ab);
      // etnte_ptr entry_ab = Props::get_non_tree_edge_entry(ab);
      // etnte_ptr entry_ba = Props::get_non_tree_edge_entry(ba);
      // etnte_ptr hdr = etnte_algo::get_header(entry_ab);
      
      // etnte_algo::erase(hdr, entry_ab);
      // etnte_algo::erase(hdr, entry_ba);

      Props::set_null_entry(ab);
      Props::set_null_entry(ba);

      // TODO: release
      // etnte_pool_.release(entry_ab);
      // etnte_pool_.release(entry_ba);
    }    


    #ifdef ET_AS_AUG_AVL
    #else
    // edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
    //   if (
    //     et_algo::node_height(et_nt::get_parent(hdr_a)) <
    //     et_algo::node_height(et_nt::get_parent(hdr_b)))
    //     for (et_ptr et_entry : range<et_nt>(hdr_a)) {
    //       if (dc_properties::is_loop_edge(et_entry)) {
    //         for (edge_t ab : graph_->edges(dc_properties::get_vertex(et_entry)))
    //           if (!dc_properties::is_null_entry(ab) && !dc_properties::is_tree_edge(ab))
    //             if (et_algo::get_header(Props::get_entry(ab->v1)) == hdr_b)
    //               return ab;
    //       }
    //     }
    //   else
    //     for (et_ptr et_entry : range<et_nt>(hdr_b)) {
    //       if (dc_properties::is_loop_edge(et_entry)) {
    //         for (edge_t ba : graph_->edges(dc_properties::get_vertex(et_entry)))
    //           if (!dc_properties::is_null_entry(ba) && !dc_properties::is_tree_edge(ba))
    //             if (et_algo::get_header(Props::get_entry(ba->v1)) == hdr_a)
    //               return Props::get_opposite(ba);
    //       }
    //     }
    //
    //
    //   return nullptr;
    // }



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

    edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
      auto [iter_a, end_a] = range<et_nt>(hdr_a);
      auto [iter_b, end_b] = range<et_nt>(hdr_b);
    
      while (true) {
        if (iter_a == end_a || iter_b == end_b)
          return nullptr;
    
        if (dc_properties::is_loop_edge(*iter_a)) {
          for (edge_t ab : graph_->edges(dc_properties::get_vertex(*iter_a)))
            if (!dc_properties::is_null_entry(ab) && !dc_properties::is_tree_edge(ab))
              if (et_algo::get_header(Props::get_entry(ab->v1)) == hdr_b)
                return ab;
        }
    
        if (dc_properties::is_loop_edge(*iter_b)) {
          for (edge_t ba : graph_->edges(dc_properties::get_vertex(*iter_b)))
            if (!dc_properties::is_null_entry(ba) && !dc_properties::is_tree_edge(ba))
              if (et_algo::get_header(Props::get_entry(ba->v1)) == hdr_a)
                return Props::get_opposite(ba);
        }
    
        ++iter_a;
        ++iter_b;
      }
    }
    #endif

    /**
     * \return true if reconnected
     */
    bool tree_edge_split_and_reconnect(edge_t ab) { // TODO: check release of entries
      auto ba = Props::get_opposite(ab);

      et_ptr et_ab = Props::get_tree_edge_entry(ab);
      et_ptr et_ba = Props::get_tree_edge_entry(ba);
      
      et_ptr et_hdr_a = et_algo::get_header(et_ab);
      // etnte_ptr etnte_hdr_a = Props::get_etnte_header(et_hdr_a);

      et_ptr et_hdr_b = et_pool_.acquire();
      // etnte_ptr etnte_hdr_b = etnte_pool_.acquire();
      et_algo::init_header(et_hdr_b);
      // etnte_algo::init_header(etnte_hdr_b);                           
      // Props::set_etnte_header(et_hdr_b, etnte_hdr_b);                        

      // etnte_context_operations<Props>::split(etnte_hdr_a, etnte_hdr_b, et_ab, et_ba);

      if (et_algo::less_than(et_ab, et_ba))
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ab, et_ba);
      else {
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ba, et_ab);
        et_algo::swap_tree(et_hdr_a, et_hdr_b);        
      }      
                                                                
      Props::set_null_entry(ab);
      Props::set_null_entry(ba);

      // bool swap = Props::less_than(etnte_hdr_b, etnte_hdr_a);

      // if (swap) {
      //   std::swap(et_hdr_a, et_hdr_b);
      //   std::swap(etnte_hdr_a, etnte_hdr_b);
      // }

      


      


      // for (et_ptr et_entry : range<et_nt>(et_hdr))
      //   if (dc_properties::is_loop_edge(et_entry))
      //     for (edge_t de : g.edges(dc_properties::get_vertex(et_entry)))
      //       if (dc_properties::is_null_entry(de)) { // non-tree edges
      //
      //
      //         nte_edges[i].push_front(de);
      //
      //
      //
      //
      //         etnte_ptr etnte_entry = etnte_pool_.acquire();
      //         dc_properties::set_non_tree_edge_entry(de, etnte_entry); 
      //         dc_properties::set_directed_edge(etnte_entry, de);
      //
      //         // etnte_algo::push_back(etnte_hdr, etnte_entry);
      //         boost::intrusive::avltree_algorithms<etnte_traits>::push_back(etnte_hdr, etnte_entry);
      //       }




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

      


      if (ab = find_replacement(et_hdr_a, et_hdr_b); ab) {                                
        // ab = Props::get_directed_edge(replace_ab_entry);
        ba = Props::get_opposite(ab);        
        // etnte_ptr replace_ba_entry = Props::get_non_tree_edge_entry(ba);

        et_ptr replace_a_entry = Props::get_ordering_vertex_entry(ab);
        et_ptr replace_b_entry = Props::get_ordering_vertex_entry(ba);                                        

        Props::set_directed_edge(et_ab, ab);
        Props::set_directed_edge(et_ba, ba);

        Props::set_tree_edge_entry(ab, et_ab);
        Props::set_tree_edge_entry(ba, et_ba);

        // etnte
        // cyclic<etnte_algo>::cut(etnte_hdr_a, etnte_context_operations<Props>::lower_bound(etnte_hdr_a, replace_a_entry));
        // auto least_b = etnte_context_operations<Props>::lower_bound(etnte_hdr_b, replace_b_entry);
        // cyclic<etnte_algo>::cut_least_dropped(etnte_hdr_b, least_b);
        // etnte_algo::join_trees(etnte_hdr_a, least_b, etnte_hdr_b);
        //
        // etnte_algo::erase(etnte_hdr_a, replace_ab_entry);
        // etnte_algo::erase(etnte_hdr_a, replace_ba_entry);    

        // et
        cyclic<et_algo>::cut(et_hdr_a, replace_a_entry);
        cyclic<et_algo>::cut(et_hdr_b, replace_b_entry);                
                  
        et_algo::join_trees(et_hdr_a, et_ab, et_hdr_b);        
        et_algo::push_back(et_hdr_a, et_ba);                       
          
        return true;
      }      

      return false;
    }


        


    void release_component(et_ptr hdr) { // TODO: check consistency
      etnte_ptr etnte_hdr = Props::get_etnte_header(hdr);

      etnte_ptr etnte_node = etnte_nt::get_left(etnte_hdr);
      while (etnte_node != etnte_hdr) {
        edge_t de = Props::get_directed_edge(etnte_node);
        Props::set_null_entry(de);
        etnte_ptr prev = etnte_node;
        etnte_node = etnte_algo::next_node(etnte_node);
        etnte_pool_.release(prev);
      }

      etnte_pool_.release(etnte_hdr);


      et_ptr et_node = et_nt::get_left(hdr);      
      while (et_node != hdr) {
        if (Props::is_loop_edge(et_node))
          Props::set_entry(Props::get_vertex(et_node), nullptr);
        else
          Props::set_null_entry(Props::get_directed_edge(et_node));

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
