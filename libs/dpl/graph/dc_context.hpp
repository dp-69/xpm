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

  template <typename Props>
  class dc_context
  {
    using et_ptr = et_traits::node_ptr;
    using et_nt = et_traits;

    using vertex_t = dc_graph::vertex_descriptor;
    using edge_t = dc_graph::edge_descriptor;

    dc_graph* g_;
    dc_properties props_; // TODO: INIT
    // std::unique_ptr<std::list<edge_t>[]> nte_edges;

    static inline constexpr auto null_edge_ = std::numeric_limits<edge_t>::max();

  public:
    // std::vector<std::size_t> saved_replacement_indices;
    // std::size_t quick_replacement_idx = 0;

    helper::smart_pool<et_traits::node> et_pool_;

    void init_with_dfs(dc_graph& g, dc_properties props) {
      g_ = &g;
      props_ = props;

      // saved_replacement_indices.reserve(g.edge_count());

      auto vertex_count = num_vertices(g);                 

      using color_t = boost::two_bit_color_type;
      auto color_uptr = std::make_unique<color_t[]>(vertex_count);

      boost::iterator_property_map color_map{
        color_uptr.get(),
        boost::identity_property_map{}
      };

      auto tree_edge_stack = std::vector<et_ptr>(vertex_count + 1);     
      euler_tour_visitor<dc_graph, dc_properties> visitor{props_, &et_pool_, tree_edge_stack.data()};

      for (std::size_t i = 0; i < g.vertex_count(); ++i) {
        if (color_map[i] != color_t::two_bit_black) {
          depth_first_visit(g, i, visitor, color_map); // tree edges

          for (et_ptr et : range<et_nt>(et_algo::get_header(props_.get_entry(i))))
            if (dc_properties::is_loop_edge(et))
              for (edge_t de : g.edges(dc_properties::get_vertex(et)))
                if (props_.is_null_entry(de)) // non-tree edges
                  props_.set_non_tree_edge(de);
        }
      }
    }

    void adjacent_edges_remove(vertex_t v, const dc_graph& g) {
      for (edge_t de : g.edges(v))
        if (!props_.is_null_entry(de) && !props_.is_tree_edge(de))
          non_tree_edge_remove(de);

      for (edge_t de : g.edges(v))
        if (!props_.is_null_entry(de))
          tree_edge_split_and_reconnect(de);

      props_.set_entry(v, nullptr);
      // TODO: release vertex?
    }

    void non_tree_edge_remove(edge_t ab) {
      edge_t ba = opposite(ab, *g_);
      props_.set_null_entry(ab);
      props_.set_null_entry(ba);
    }    


    edge_t find_replacement(et_ptr hdr_a, et_ptr hdr_b) {
      if (
        et_algo::node_height(et_nt::get_parent(hdr_a)) <
        et_algo::node_height(et_nt::get_parent(hdr_b)))
        for (et_ptr et_entry : range<et_nt>(hdr_a)) {
          if (dc_properties::is_loop_edge(et_entry)) {
            for (edge_t ab : g_->edges(dc_properties::get_vertex(et_entry)))
              if (!props_.is_null_entry(ab) && !props_.is_tree_edge(ab))
                if (et_algo::get_header(props_.get_entry(target(ab, *g_))) == hdr_b)
                  return ab;
          }
        }
      else
        for (et_ptr et_entry : range<et_nt>(hdr_b)) {
          if (dc_properties::is_loop_edge(et_entry)) {
            for (edge_t ba : g_->edges(dc_properties::get_vertex(et_entry)))
              if (!props_.is_null_entry(ba) && !props_.is_tree_edge(ba))
                if (et_algo::get_header(props_.get_entry(target(ba, *g_))) == hdr_a)
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

      et_ptr et_ab = props_.get_tree_edge_entry(ab);
      et_ptr et_ba = props_.get_tree_edge_entry(ba);
      
      et_ptr et_hdr_a = et_algo::get_header(et_ab);
      et_ptr et_hdr_b = et_pool_.acquire();
      et_algo::init_header(et_hdr_b);

      if (et_algo::less_than(et_ab, et_ba))
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ab, et_ba);
      else {
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ba, et_ab);
        et_algo::swap_tree(et_hdr_a, et_hdr_b);        
      }      
                                                                
      props_.set_null_entry(ab);
      props_.set_null_entry(ba);





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

        et_ptr replace_a_entry = props_.get_ordering_vertex_entry(ab);
        et_ptr replace_b_entry = props_.get_ordering_vertex_entry(ba);                                        

        Props::set_directed_edge(et_ab, ab);
        Props::set_directed_edge(et_ba, ba);

        props_.set_tree_edge_entry(ab, et_ab);
        props_.set_tree_edge_entry(ba, et_ba);

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
