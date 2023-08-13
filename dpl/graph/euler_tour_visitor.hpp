#pragma once


#include "et_etnte_defs.hpp"
#include "smart_pool.hpp"
#include <stack>

//#include "HW/pore_network_modelling/pnm_static.hpp"

namespace HW { namespace dynamic_connectivity
{

  class et_context
  {
    using et_node_ptr = et_traits::node_ptr;
    using et_const_node_ptr = et_traits::const_node_ptr;

    using etnte_node = dpl::graph::aug_avl_node;

    using mask = dpl::graph::mask;

  public:
    static etnte_node* get_non_tree_edge_header(et_const_node_ptr n) {
      return mask::get_ptr<etnte_node>(n->tag);
    }

    static void set_non_tree_edge_header(et_node_ptr n, const etnte_node* etnte_header) {
      mask::set_ptr_balance(n->tag, etnte_header, mask::balance);
    }

    static directed_edge* get_directed_edge(et_const_node_ptr n) {
      return mask::get_ptr<directed_edge>(n->tag);
    }

    static vertex* get_vertex(et_const_node_ptr n) {
      return mask::get_ptr<vertex>(n->tag);
    }
    
    static void set_directed_edge(et_node_ptr n, const directed_edge* de) {
      mask::set_ptr(n->tag, de);
    }
    
    static void set_vertex(et_node_ptr n, const vertex* v) {
      mask::set_ptr_bit(n->tag, v);
    }                     

    static bool is_loop_edge(et_const_node_ptr n) { // that is 'vertex'
      return mask::get_bit(n->tag);
    }
  };





  struct etnte_context
  {
    static directed_edge* get_directed_edge(etnte_traits::const_node_ptr n) {
      return dpl::graph::mask::get_ptr<directed_edge>(n->tag);
    }

    static void set_directed_edge(etnte_traits::node_ptr n, const directed_edge* de) {
      dpl::graph::mask::set_ptr(n->tag, de);
    }

    // ordering of non-tree edges
    static et_traits::node_ptr get_vertex_entry(const directed_edge* de) {
      // sorted by a pointing-in vertex, the pointing-out works as well        
      return de->opposite->v1->et_entry_;   
    }
    
    static et_traits::node_ptr get_vertex_entry(etnte_traits::const_node_ptr n) {
      return get_vertex_entry(get_directed_edge(n));      
    }
  };


































//  struct component_description
//  {
//    et_node_ptr header,
//    size_t vertexCount,
//    size_t
//  };

//  using two_phase_flow::phase_tag;

//  class isolated_visitor : public boost::empty_dfs_visitor<decremental_connectivity_graph>
//  {
//    vector<phase_tag>& _edgeConductance;
//    const graph_storage& _storage;
//
//  public:
//    isolated_visitor(vector<phase_tag>& edgeConductance, const graph_storage& storage)
//      : _edgeConductance(edgeConductance),
//        _storage(storage) {}
//
//    void forward_or_cross_edge(directed_edge_ptr e, const decremental_connectivity_graph&) {
//      _edgeConductance[_storage.edge_idx(e)] = phase_tag::water;
//    }
//
//    void tree_edge(directed_edge_ptr e, const decremental_connectivity_graph&) {
//      _edgeConductance[_storage.edge_idx(e)] = phase_tag::water;
//    }
//  };

  

  
  class euler_tour_visitor : public boost::empty_dfs_visitor<dynamic_connectivity_graph>
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;



    et_node_ptr _etHeader;
    etnte_node_ptr _etnteHeader;
//    et_node_ptr _etPtr;
//    etnte_node_ptr _etntePtr;

    smart_pool<et_traits::node>& _etPool;
    smart_pool<etnte_traits::node>& _etntePool;

    et_node_ptr* _treeEdgeStackEmpty;
    et_node_ptr* _treeEdgeStackTop;


    std::vector<et_node_ptr>& _components;        

    et_node_ptr acquire_zero_level_et_entry() {
//      auto entry = _etPtr++;
      auto entry = _etPool.acquire();
//      et_traits::set_next_level_entry(entry, nullptr);                
      
//      et_traits::augment_identity(entry);
      return entry;
    }    
    
  public:
    

    euler_tour_visitor(
      smart_pool<et_traits::node>& etPool, smart_pool<etnte_traits::node>& etntePool, et_node_ptr* treeEdgeStack, std::vector<et_node_ptr>& components)
      : _etPool(etPool),
        _etntePool(etntePool),
        _treeEdgeStackEmpty(treeEdgeStack),
        _treeEdgeStackTop(treeEdgeStack),
        _components(components)    
    {}
    
    



    // non-tree edge    
    void forward_or_cross_edge(directed_edge_ptr de, const dynamic_connectivity_graph&) {
      auto entry = _etntePool.acquire();
      auto entryOpp = _etntePool.acquire();
      
      auto deOpp = de->opposite;

      directed_edge::set_non_tree_edge_entry(de, entry); 
      directed_edge::set_non_tree_edge_entry(deOpp, entryOpp);
      
      etnte_context::set_directed_edge(entry, de);
      etnte_context::set_directed_edge(entryOpp, deOpp);
            
      etnte_context_operations<etnte_context>::insert(_etnteHeader, entry);
      etnte_context_operations<etnte_context>::insert(_etnteHeader, entryOpp);                    
    }

    void start_vertex(vertex_ptr, const dynamic_connectivity_graph&) {            
      _etHeader = _etPool.acquire();
      et_algo::init_header(_etHeader);      
        
      _etnteHeader = _etntePool.acquire();
      etnte_algo::init_header(_etnteHeader);      
      et_context::set_non_tree_edge_header(_etHeader, _etnteHeader);             

      _components.push_back(_etHeader);
    }    

    void discover_vertex(vertex_ptr v, const dynamic_connectivity_graph&) {     
      auto entry = acquire_zero_level_et_entry();        
      
      et_context::set_vertex(entry, v);
      v->et_entry_ = entry;
//      vertex_color_property_map::set_et_entry(v, entry);

//      et_algo_not_augmented::push_back(_etHeader, entry);            
      et_algo::push_back(_etHeader, entry);            
    }

    void tree_edge(directed_edge_ptr e, const dynamic_connectivity_graph&) {
      auto entry = acquire_zero_level_et_entry();

      et_context::set_directed_edge(entry, e);
      directed_edge::set_tree_edge_entry(e, entry);

      et_algo::push_back(_etHeader, entry);      

      *++_treeEdgeStackTop = entry;
    }

    void finish_vertex(vertex_ptr v, const dynamic_connectivity_graph&) {
      if (_treeEdgeStackTop != _treeEdgeStackEmpty) {
        auto entry = acquire_zero_level_et_entry();
        
        auto top = *_treeEdgeStackTop--;

        auto top_de_opposite = et_context::get_directed_edge(top)->opposite;
        et_context::set_directed_edge(entry, top_de_opposite);        
        directed_edge::set_tree_edge_entry(top_de_opposite, entry);

//        et_algo_not_augmented::push_back(_etHeader, entry);        
        et_algo::push_back(_etHeader, entry);        
      }
//      else
//        et_algo::refresh_size(_etHeader);
    }    
  }; 
  
  template<class DFSVisitor>
  static void execute_dfs(const dynamic_connectivity_graph& g, vertex_ptr rootVertex, DFSVisitor visitor) {
    boost::graph_traits<dynamic_connectivity_graph>::vertex_iterator vi, vi_start, vi_end;
    boost::tie(vi_start, vi_end) = vertices(g);
            
    for (vi = vi_start; vi != vi_end; ++vi)
      vertex_color_property_map::compress_and_init(*vi);

    depth_first_visit(g, rootVertex, visitor, vertex_color_property_map());

//      depth_first_search(g,
//        visitor(euler_tour_visitor(etPool, etntePool, treeEdgeStack.data(), components))
//          .color_map(vertex_color_property_map())/*.root_vertex(rootVertex)*/);

    for (vi = vi_start; vi != vi_end; ++vi)
      vertex_color_property_map::decompress_and_finish(*vi);
  }


















  


  static bool execute_connectivity_dfs(dynamic_connectivity_graph& g, vertex_ptr searchVertex, vertex_ptr outletVertex) {    
//    if (searchVertex->trapped)
//      return false;
    

    boost::graph_traits<dynamic_connectivity_graph>::vertex_iterator vi, vi_start, vi_end;
    boost::tie(vi_start, vi_end) = vertices(g);
            
    for (vi = vi_start; vi != vi_end; ++vi)
      vi->visited = false;

    
    std::stack<vertex_ptr> vertexStack;    
    vertexStack.push(searchVertex);

    while (!vertexStack.empty()) {
      const auto u = vertexStack.top();
      vertexStack.pop();
      if (!u->visited) {
        u->visited = true;

        // Discovery        
        if (u == outletVertex)
           return true;
        

        for (const auto& de : out_edges_range(u))
          if (!de->v1->visited)
            vertexStack.push(de->v1);
      }
    }


    

    for (vi = vi_start; vi != vi_end; ) {
      auto prev_v = *vi;
      ++vi;
      if (prev_v->visited) {        
//        clear_out_edges(prev_v, g);
        remove_vertex(prev_v, g);      
      }
    }



    return false;


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
  }

}}
