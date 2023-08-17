#pragma once


#include "et_etnte_defs.hpp"
#include <stack>

namespace HW::dynamic_connectivity
{
  class etnte_context
  {
    using et_node_ptr = et_traits::node_ptr;
    using et_const_node_ptr = et_traits::const_node_ptr;

    using etnte_node = etnte_traits::node;
    using etnte_node_ptr = etnte_traits::node_ptr;
    using etnte_const_node_ptr = etnte_traits::const_node_ptr;

    using mask = dpl::graph::mask;

  public:
    static etnte_node_ptr get_non_tree_edge_header(et_const_node_ptr n) {
      return mask::get_ptr<etnte_node>(n->tag);
    }

    static void set_non_tree_edge_header(et_node_ptr n, etnte_const_node_ptr etnte_header) {
      mask::set_ptr_balance(n->tag, etnte_header, mask::balance);
    }

    static directed_edge* get_directed_edge(et_const_node_ptr n) {
      return mask::get_ptr<directed_edge>(n->tag);
    }

    static vertex_ptr get_vertex(et_const_node_ptr n) {
      return mask::get_ptr<vertex>(n->tag);
    }
    
    static void set_directed_edge(et_node_ptr n, const directed_edge_ptr de) {
      mask::set_ptr(n->tag, de);
    }
    
    static void set_vertex(et_node_ptr n, const vertex_ptr v) {
      mask::set_ptr_bit(n->tag, v);
    }                     

    static bool is_loop_edge(et_const_node_ptr n) { // that is 'vertex'
      return mask::get_bit(n->tag);
    }

    // ---------------



    static directed_edge_ptr get_directed_edge(etnte_const_node_ptr n) {
      return mask::get_ptr<directed_edge>(n->tag);
    }

    static void set_directed_edge(etnte_node_ptr n, const directed_edge_ptr de) {
      mask::set_ptr(n->tag, de);
    }

    // ordering of non-tree edges
    static et_node_ptr get_ordering_vertex_entry(const directed_edge_ptr de) {
      // sorted by a pointing-in vertex, the pointing-out works as well        
      return de->opposite->v1->et_entry_;   
    }
    
    static et_node_ptr get_ordering_vertex_entry(etnte_const_node_ptr n) {
      return get_ordering_vertex_entry(get_directed_edge(n));      
    }

    // ---------------

    static et_node_ptr get_entry(vertex_ptr v) {
      return v->et_entry_;
    }

    static void set_entry(vertex_ptr v, et_node_ptr et) {
      v->et_entry_ = et;
    }


    // ---------------

    typedef tagged_pointer_as_size_t<bool, 1> compression;

    static bool is_tree_edge(const directed_edge_ptr& x) {
      return compression::get_bits(x->entry_type_);
    }

    static et_node_ptr get_tree_edge_entry(const directed_edge_ptr& x) {
      return compression::get_pointer<et_node_ptr>(x->entry_type_);      
    }

    static void set_tree_edge_entry(const directed_edge_ptr& x, const et_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, true);      
    }

    static etnte_node_ptr get_non_tree_edge_entry(const directed_edge_ptr& x) {
      return compression::get_pointer<etnte_node_ptr>(x->entry_type_);      
    }

    static void set_non_tree_edge_entry(const directed_edge_ptr& x, const etnte_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, false);      
    }    
  };

































  

  template <typename Context>
  class euler_tour_visitor : public boost::default_dfs_visitor
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    et_node_ptr et_hdr_;
    etnte_node_ptr etnte_hdr_;

    dpl::graph::smart_pool<et_traits::node>* et_pool_;
    dpl::graph::smart_pool<etnte_traits::node>* etnte_pool_;

    et_node_ptr* tree_edge_stack_empty_;
    et_node_ptr* tree_edge_stack_top_;

    // std::vector<et_node_ptr>& components_;        

    
  public:
    euler_tour_visitor(
      dpl::graph::smart_pool<et_traits::node>* et_pool
    , dpl::graph::smart_pool<etnte_traits::node>* etnte_pool
    , et_node_ptr* tree_edge_stack
    )
      : et_pool_(et_pool)
      , etnte_pool_(etnte_pool)
      , tree_edge_stack_empty_(tree_edge_stack)
      , tree_edge_stack_top_(tree_edge_stack)
    {}


    void start_vertex(vertex_ptr, const dc_graph&) {            
      et_hdr_ = et_pool_->acquire();
      et_algo::init_header(et_hdr_);      
      etnte_hdr_ = etnte_pool_->acquire();
      etnte_algo::init_header(etnte_hdr_);      
      Context::set_non_tree_edge_header(et_hdr_, etnte_hdr_);             
    }

    void discover_vertex(vertex_ptr v, const dc_graph&) {     
      et_node_ptr entry = et_pool_->acquire();        
      Context::set_vertex(entry, v);
      Context::set_entry(v, entry);
      et_algo::push_back(et_hdr_, entry);            
    }

    void finish_vertex(vertex_ptr v, const dc_graph&) {
      if (tree_edge_stack_top_ != tree_edge_stack_empty_) {
        et_node_ptr entry = et_pool_->acquire();
        et_node_ptr top = *tree_edge_stack_top_--;
        directed_edge_ptr top_de_opposite = Context::get_directed_edge(top)->opposite;
        Context::set_directed_edge(entry, top_de_opposite);        
        Context::set_tree_edge_entry(top_de_opposite, entry);
        et_algo::push_back(et_hdr_, entry);        
      }
    }    

    void tree_edge(directed_edge_ptr e, const dc_graph&) {
      et_node_ptr entry = et_pool_->acquire();
      Context::set_directed_edge(entry, e);
      Context::set_tree_edge_entry(e, entry);
      et_algo::push_back(et_hdr_, entry);      
      *++tree_edge_stack_top_ = entry;
    }

    // non-tree edge    
    void forward_or_cross_edge(directed_edge_ptr de, const dc_graph&) {
      etnte_node_ptr entry = etnte_pool_->acquire();
      etnte_node_ptr entry_opp = etnte_pool_->acquire();
      
      directed_edge_ptr de_opp = de->opposite;

      Context::set_non_tree_edge_entry(de, entry); 
      Context::set_non_tree_edge_entry(de_opp, entry_opp);
      
      Context::set_directed_edge(entry, de);
      Context::set_directed_edge(entry_opp, de_opp);
            
      etnte_context_operations<Context>::insert(etnte_hdr_, entry);
      etnte_context_operations<Context>::insert(etnte_hdr_, entry_opp);                    
    }
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

  template <typename Context>
  static void execute_dfs(const dc_graph& g, euler_tour_visitor<Context> visitor) {
    vertex_iterator vi, vi_start, vi_end;
    boost::tie(vi_start, vi_end) = vertices(g);
            
    for (vi = vi_start; vi != vi_end; ++vi)
      vertex_color_property_map::compress_and_init(*vi);

    for (vi = vi_start; vi != vi_end; ++vi) {
      if (vertex_color_property_map::get_color(*vi) != boost::two_bit_color_type::two_bit_black)
        depth_first_visit(g, *vi, visitor, vertex_color_property_map());
    }

    for (vi = vi_start; vi != vi_end; ++vi)
      vertex_color_property_map::decompress_and_finish(*vi);
  }
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