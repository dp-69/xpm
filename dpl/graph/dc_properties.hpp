#pragma once

#include "et_etnte_defs.hpp"
#include "dc_graph.hpp"

namespace dpl::graph
{
  class dc_properties
  {
    using et_ptr = et_traits::node_ptr;
    using et_cptr = et_traits::const_node_ptr;

    using etnte = etnte_traits::node;
    using etnte_ptr = etnte_traits::node_ptr;
    using etnte_cptr = etnte_traits::const_node_ptr;

    const dc_graph* g_;

    using vertex_desc = dc_graph::vertex_descriptor;

  public:
    dc_properties() = default;

    explicit dc_properties(const dc_graph& graph)
      : g_(&graph) {}

    static directed_edge* get_directed_edge(et_cptr n) {
      return mask_bit_balance::get_ptr<directed_edge*>(n->tag);
    }

    static void set_directed_edge(et_ptr n, const directed_edge* de) {
      mask_bit_balance::set_ptr(n->tag, de);
    }

    static vertex_desc get_vertex(et_cptr n) {
      return mask_bit_balance::get_value<vertex_desc>(n->tag);
    }

    static void set_vertex(et_ptr n, vertex_desc v) {
      mask_bit_balance::set_value_bit(n->tag, v);
    }                     

    /**
     * \brief checks if n is a vertex-type et entry
     */
    static bool is_loop_edge(et_cptr n) {
      return mask_bit_balance::get_bit(n->tag);
    }

    // ---------------

    /**
     * \brief
     *   ordering of non-tree edges
     *   sorted by a pointing-in vertex, the pointing-out works as well
     */
    et_ptr get_ordering_vertex_entry(const directed_edge* de) const {
      return get_entry(target(opposite(de, *g_), *g_));
    }

    // ---------------

    et_ptr get_entry(vertex_desc v) const {
      return static_cast<et_ptr>(g_->vertex_entry(v));
    }

    void set_entry(vertex_desc v, et_ptr et) {
      g_->vertex_entry(v) = et;
    }

    // ---------------

    bool is_tree_edge(directed_edge* x) const {
      return !mask_bit::get_bit(g_->directed_edge_entry(x));
    }

    et_ptr get_tree_edge_entry(directed_edge* x) const { return mask_bit::get_ptr<et_ptr>(g_->directed_edge_entry(x)); }
    void set_tree_edge_entry(directed_edge* x, et_ptr y) { mask_bit::set_ptr(g_->directed_edge_entry(x), y); }

    // static etnte_ptr get_non_tree_edge_entry(const directed_edge* x) { return mask_bit::get_ptr<etnte_ptr>(x->entry_); }
    void set_non_tree_edge(directed_edge* x) { return mask_bit::set_bit(g_->directed_edge_entry(x)); }

    bool is_null_entry(directed_edge* x) const { return g_->directed_edge_entry(x) == 0; }
    void set_null_entry(directed_edge* x) { g_->directed_edge_entry(x) = 0; }

    // static bool less_than(etnte_cptr hdr_l, etnte_cptr hdr_r) {
    //   etnte_ptr root_l = etnte_traits::get_parent(hdr_l);
    //   etnte_ptr root_r = etnte_traits::get_parent(hdr_r);
    //
    //   return
    //     #ifdef ETNTE_AS_AVL_ONLY
    //       (root_l ? etnte_algo::node_height(root_l) : 0) < 
    //       (root_r ? etnte_algo::node_height(root_r) : 0);
    //     #else
    //       (root_l ? etnte_traits::get_size(root_l) : 0) < 
    //       (root_r ? etnte_traits::get_size(root_r) : 0);
    //     #endif
    // }

    // static bool less_than(etnte_cptr hdr_l, etnte_cptr hdr_r) {
    //   etnte_ptr root_l = etnte_traits::get_parent(hdr_l);
    //   etnte_ptr root_r = etnte_traits::get_parent(hdr_r);
    //
    //   return
    //     #ifdef ETNTE_AS_AVL_ONLY
    //       (root_l ? etnte_algo::node_height(root_l) : 0) < 
    //       (root_r ? etnte_algo::node_height(root_r) : 0);
    //     #else
    //       (root_l ? etnte_traits::get_size(root_l) : 0) < 
    //       (root_r ? etnte_traits::get_size(root_r) : 0);
    //     #endif
    // }
  };












  

  
  
  // ReSharper disable once CppPassValueParameterByConstReference
  // static void execute_dfs(const dynamic_connectivity_graph& g, vertex_ptr root, euler_tour_visitor visitor) {
  //   boost::graph_traits<dynamic_connectivity_graph>::vertex_iterator vi, vi_start, vi_end;
  //   boost::tie(vi_start, vi_end) = vertices(g);
  //           
  //   for (vi = vi_start; vi != vi_end; ++vi)
  //     vertex_color_property_map::compress_and_init(*vi);
  //
  //   depth_first_visit(g, root, visitor, vertex_color_property_map());
  //
  //   for (vi = vi_start; vi != vi_end; ++vi)
  //     vertex_color_property_map::decompress_and_finish(*vi);
  // }

  
}





// static bool execute_connectivity_dfs(dynamic_connectivity_graph& g, vertex_ptr searchVertex, vertex_ptr outletVertex) {    
// //    if (searchVertex->trapped)
// //      return false;
//     
//
//     boost::graph_traits<dynamic_connectivity_graph>::vertex_iterator vi, vi_start, vi_end;
//     boost::tie(vi_start, vi_end) = vertices(g);
//             
//     for (vi = vi_start; vi != vi_end; ++vi)
//       vi->visited = false;
//
//     
//     std::stack<vertex_ptr> vertexStack;    
//     vertexStack.push(searchVertex);
//
//     while (!vertexStack.empty()) {
//       const auto u = vertexStack.top();
//       vertexStack.pop();
//       if (!u->visited) {
//         u->visited = true;
//
//         // Discovery        
//         if (u == outletVertex)
//            return true;
//         
//
//         for (const auto& de : out_edges_range(u))
//           if (!de->v1->visited)
//             vertexStack.push(de->v1);
//       }
//     }
//
//
//     
//
//     for (vi = vi_start; vi != vi_end; ) {
//       auto prev_v = *vi;
//       ++vi;
//       if (prev_v->visited) {        
// //        clear_out_edges(prev_v, g);
//         remove_vertex(prev_v, g);      
//       }
//     }
//
//
//
//     return false;
//   }



























//    boost::depth_first_visit(g, searchVertex, connectivity_visitor(outletVertex), vertex_color_property_map());
//
//
//
//
//    auto found = false;
//
//    auto vis = connectivity_visitor(outletVertex);
//    auto u = searchVertex;
//    auto color = vertex_color_property_map();
//
//    vis.start_vertex(u, g);
//    typedef dynamic_connectivity_graph IncidenceGraph;
//  
//      typedef boost::graph_traits<IncidenceGraph>::vertex_descriptor Vertex;
//      typedef boost::graph_traits<IncidenceGraph>::edge_descriptor Edge;
// 
//      typedef boost::property_traits<vertex_color_property_map>::value_type ColorValue;
//  
//      typedef boost::color_traits<ColorValue> Color;
//      typedef boost::graph_traits<IncidenceGraph>::out_edge_iterator Iter;
//      typedef pair<Vertex, pair<boost::optional<Edge>, pair<Iter, Iter> > > VertexInfo;
//
//      boost::optional<Edge> src_e;
//      Iter ei, ei_end;
//      vector<VertexInfo> stack;
//
//      // Possible optimization for vector
//      //stack.reserve(num_vertices(g));
//
//      put(color, u, Color::gray());
//      vis.discover_vertex(u, g);
//      if (u == outletVertex) {
//        found = true;
//        goto end;
//      }
//
//
//      boost::tie(ei, ei_end) = out_edges(u, g);
//      if (/*func(u, g)*/ false) {
//          // If this vertex terminates the search, we push empty range
//          stack.push_back(std::make_pair(u, std::make_pair(boost::optional<Edge>(), std::make_pair(ei_end, ei_end))));
//      } else {
//          stack.push_back(std::make_pair(u, std::make_pair(boost::optional<Edge>(), std::make_pair(ei, ei_end))));
//      }
//      while (!stack.empty()) {
//        VertexInfo& back = stack.back();
//        u = back.first;
//        src_e = back.second.first;
//        boost::tie(ei, ei_end) = back.second.second;
//        stack.pop_back();
//	// finish_edge has to be called here, not after the
//	// loop. Think of the pop as the return from a recursive call.
////        if (src_e) {
////	  boost::call_finish_edge(vis, src_e.get(), g);
////	}
//        while (ei != ei_end) {
//          Vertex v = target(*ei, g);
//          vis.examine_edge(*ei, g);
//          ColorValue v_color = get(color, v);
//          if (v_color == Color::white()) {
//            vis.tree_edge(*ei, g);
//            src_e = *ei;
//            stack.push_back(std::make_pair(u, std::make_pair(src_e, std::make_pair(++ei, ei_end))));
//            u = v;
//            put(color, u, Color::gray());
//            vis.discover_vertex(u, g);
//            if (u == outletVertex) {
//              found = true;
//              goto end;
//            }
//
//            boost::tie(ei, ei_end) = out_edges(u, g);
//            if (/*func(u, g)*/ false) {
//                ei = ei_end;
//            }
//          } else {
//            if (v_color == Color::gray()) {
//              vis.back_edge(*ei, g);
//            } else {
//              vis.forward_or_cross_edge(*ei, g);
//            }
////            call_finish_edge(vis, *ei, g);
//            ++ei;
//          }
//        }
//        put(color, u, Color::black());
//        vis.finish_vertex(u, g);
//      }
//
//
//
//      end:
//
//      for (vi = vi_start; vi != vi_end; ++vi) {
//        if (!found && get(color, *vi) != Color::white())
//          vi->trapped = true;
//
//        vertex_color_property_map::decompress_and_finish(*vi);
//      }
//
//
//
//      return found;




//     depth_first_visit(g, searchVertex, connectivity_visitor(outletVertex), vertex_color_property_map());
//    auto found = false;
//
//    try {
//      depth_first_visit(g, rootVertex, connectivity_visitor(outletVertex), vertex_color_property_map());
//    }
//    catch(int) {
//      found = true;
//    }




//      depth_first_search(g,
//        visitor(euler_tour_visitor(etPool, etntePool, treeEdgeStack.data(), components))
//          .color_map(vertex_color_property_map())/*.root_vertex(rootVertex)*/);

//    auto connected = outletVertex->connectState == 1;
//
//    for (vi = vi_start; vi != vi_end; ++vi) {
//      vertex_color_property_map::decompress_and_finish(*vi);
//      if (vi->connectState == 1) {
//        if (connected)
//          vi->connectState = 0;
//        else {
//          vi->connectState = 2; // trapped
//        }
//      }
//      
//
//    }
//
//    return connected;





//    return found;