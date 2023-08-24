#pragma once

#include "dc_properties.hpp"

#include <boost/graph/two_bit_color_map.hpp>
#include <boost/property_map/function_property_map.hpp>

namespace dpl::graph
{      
  struct component_visitor_empty
  {
    void vertex(vertex*) {}
    void directed_edge(directed_edge*) {}
  };
    
  struct component_visitor_count
  {
    size_t _vertices = 0;
    size_t _directedEdges = 0;

    void vertex(vertex*) { ++_vertices; }
    void directed_edge(directed_edge*) { ++_directedEdges; }

    size_t vertices() { return _vertices; }
    size_t edges() { return _directedEdges/2; }
  };


  struct dc_stat
  {
    size_t nteCountOther;
    size_t nteCountActive;
    size_t nteReconnectingFirst;
  };

  template <typename Props>
  class dc_context
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    using et_nt = et_traits;
    using etnte_nt = etnte_traits;

  public:    
    helper::smart_pool<et_traits::node> et_pool_;
    helper::smart_pool<etnte_traits::node> etnte_pool_;

    // std::vector<et_node_ptr> initial_components_;
   
    // template<class Visitor>
    // void init(const dynamic_connectivity_graph& g, vertex_ptr rootVertex, Visitor v) {            
    //   auto vertexCount = num_vertices(g);                 
    //   auto treeEdgeStack = std::vector<et_node_ptr>(vertexCount + 1);     
    //       
    //   execute_dfs(g, rootVertex,              
    //     merge_visitors(euler_tour_visitor(&etPool_, &etntePool_, treeEdgeStack.data()/*, initial_components_*/), v));
    // }
    
    // void init(const dynamic_connectivity_graph& g, vertex_ptr rootVertex) {           
    //   auto vertexCount = num_vertices(g);                 
    //   auto treeEdgeStack = std::vector<et_node_ptr>(vertexCount + 1);     
    //       
    //   execute_dfs(g, rootVertex, euler_tour_visitor(&etPool_, &etntePool_, treeEdgeStack.data()/*, initial_components_*/));
    // }

    void init(dc_graph& g, Props& c) {           
      auto vertex_count = num_vertices(g);                 

      using color_t = boost::two_bit_color_type;
      auto color_uptr = std::make_unique<color_t[]>(vertex_count);
      boost::iterator_property_map color_map{
        color_uptr.get(),
        boost::make_function_property_map<vertex*>(
          [&c](vertex* v) { return c.get_idx(v); })
      };

      auto tree_edge_stack = std::vector<et_node_ptr>(vertex_count + 1);     
      euler_tour_visitor<dc_graph, Props> visitor{&et_pool_, &etnte_pool_, tree_edge_stack.data()};

      for (vertex* v : range(g))
        if (color_map[v] != color_t::two_bit_black)
          boost::depth_first_visit(g, v, visitor, color_map);
    }
    

    void non_tree_edge_remove(directed_edge* ab) {
      auto ba = Props::get_opposite(ab);
      auto etnteAB = Props::get_non_tree_edge_entry(ab);
      auto etnteBA = Props::get_non_tree_edge_entry(ba);
      auto etnteHeader = etnte_algo::get_header(etnteAB);
      
      etnte_algo::erase(etnteHeader, etnteAB);
      etnte_algo::erase(etnteHeader, etnteBA);

      Props::set_null_entry(ab);
      Props::set_null_entry(ba);

      etnte_pool_.release(etnteAB);
      etnte_pool_.release(etnteBA);
    }    

    // Returns true if reconnected.
    bool tree_edge_split_and_reconnect(directed_edge* ab) {
      auto ba = Props::get_opposite(ab);

      et_node_ptr et_ab = Props::get_tree_edge_entry(ab);
      et_node_ptr et_ba = Props::get_tree_edge_entry(ba);
      
      et_node_ptr et_hdr_a = et_algo::get_header(et_ab);
      etnte_node_ptr etnte_hdr_a = Props::get_etnte_header(et_hdr_a);

      et_node_ptr et_hdr_b = et_pool_.acquire();
      etnte_node_ptr etnte_hdr_b = etnte_pool_.acquire();
      et_algo::init_header(et_hdr_b);
      etnte_algo::init_header(etnte_hdr_b);                           
      Props::set_etnte_header(et_hdr_b, etnte_hdr_b);                        

      etnte_context_operations<Props>::split(etnte_hdr_a, etnte_hdr_b, et_ab, et_ba);

      if (et_algo::less_than(et_ab, et_ba))
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ab, et_ba);
      else {
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ba, et_ab);
        et_algo::swap_tree(et_hdr_a, et_hdr_b);        
      }      
      
                                                                
      Props::set_null_entry(ab);
      Props::set_null_entry(ba);      


      {
        etnte_node_ptr root_a = etnte_nt::get_parent(etnte_hdr_a);
        etnte_node_ptr root_b = etnte_nt::get_parent(etnte_hdr_b);

        auto size_a = root_a ? etnte_nt::get_size(root_a) : 0;
        auto size_b = root_b ? etnte_nt::get_size(root_b) : 0;

        if (size_b < size_a) {
          std::swap(et_hdr_a, et_hdr_b);
          std::swap(etnte_hdr_a, etnte_hdr_b);
        }
      }

      etnte_node_ptr replacement_ab = nullptr;        
      
      for (etnte_node_ptr etnte : range<etnte_traits>(etnte_hdr_a))
        if (etnte_algo::get_header(
          Props::get_non_tree_edge_entry(
            Props::get_opposite(
              Props::get_directed_edge(etnte)))) == etnte_hdr_b) {
          replacement_ab = etnte;
          break;
        }
     
      if (replacement_ab) {                                
        ab = Props::get_directed_edge(replacement_ab);
        ba = Props::get_opposite(ab);        
        etnte_node_ptr replacement_ba = Props::get_non_tree_edge_entry(ba);

        et_node_ptr et_rec_a = Props::get_ordering_vertex_entry(ab);
        et_node_ptr et_rec_b = Props::get_ordering_vertex_entry(ba);                                        

        Props::set_directed_edge(et_ab, ab);
        Props::set_directed_edge(et_ba, ba);

        Props::set_tree_edge_entry(ab, et_ab);
        Props::set_tree_edge_entry(ba, et_ba);

        // etnte
        cyclic<etnte_algo>::principal_cut(etnte_hdr_a, etnte_context_operations<Props>::lower_bound(etnte_hdr_a, et_rec_a));
        auto leastB = etnte_context_operations<Props>::lower_bound(etnte_hdr_b, et_rec_b);
        cyclic<etnte_algo>::principal_cut_least_dropped(etnte_hdr_b, leastB);
        etnte_algo::join_trees(etnte_hdr_a, leastB, etnte_hdr_b);

        etnte_algo::erase(etnte_hdr_a, replacement_ab);
        etnte_algo::erase(etnte_hdr_a, replacement_ba);    

        // et
        cyclic<et_algo>::principal_cut(et_hdr_a, et_rec_a);
        cyclic<et_algo>::principal_cut(et_hdr_b, et_rec_b);                
                  
        et_algo::join_trees(et_hdr_a, et_ab, et_hdr_b);        
        et_algo::push_back(et_hdr_a, et_ba);                       
          
        return true;
      }      
     

      return false;
    }


        


    template<class Visitor = component_visitor_empty>
    void remove_and_release_component(et_node_ptr header, Visitor& visitor = Visitor()) {            
      auto etnteHeader = Props::get_etnte_header(header);

      auto etnteNode = etnte_nt::get_left(etnteHeader);
      while (etnteNode != etnteHeader) {
        auto de = Props::get_directed_edge(etnteNode);
        Props::set_null_entry(de);
        auto prev = etnteNode;
        etnteNode = etnte_algo::next_node(etnteNode);
        etnte_pool_.release(prev);
      }

      etnte_pool_.release(etnteHeader);


      auto etNode = et_nt::get_left(header);      
      while (etNode != header) {
        if (Props::is_loop_edge(etNode))
          Props::set_entry(Props::get_vertex(etNode), nullptr);
        else
          Props::set_null_entry(Props::get_directed_edge(etNode));

        etNode = et_algo::next_node(etNode);
      }

      etNode = et_nt::get_left(header);
      while (etNode != header) {
        if (Props::is_loop_edge(etNode))
          visitor.vertex(Props::get_vertex(etNode));        

        auto prev = etNode;
        etNode = et_algo::next_node(etNode);
        et_pool_.release(prev);
      }

      et_pool_.release(header);
    }


    
  };  
}
