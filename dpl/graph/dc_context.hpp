#pragma once

#include "dc_graph.hpp"
#include "dc_et_context.hpp"

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

  template <typename Context>
  class dc_context
  {
    static constexpr auto stack_capacity = 256;

    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    using et_nt = et_traits;
    using etnte_nt = etnte_traits;

  public:    

    smart_pool<et_traits::node> etPool_;
    smart_pool<etnte_traits::node> etntePool_;

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

    void init(dc_graph& g, Context& c) {           
      auto vertex_count = num_vertices(g);                 

      

      using color_t = boost::two_bit_color_type;
      auto color_uptr = std::make_unique<color_t[]>(vertex_count);
      boost::iterator_property_map color_map{
        color_uptr.get(),
        boost::make_function_property_map<vertex*>(
          [&c](vertex* v) { return c.get_idx(v); })
      };

      auto tree_edge_stack = std::vector<et_node_ptr>(vertex_count + 1);     
      euler_tour_visitor<dc_graph, Context> visitor{&etPool_, &etntePool_, tree_edge_stack.data()};



      for (vertex& v : range(g))
        if (color_map[&v] != color_t::two_bit_black)
          boost::depth_first_visit(g, &v, visitor, color_map);
    }

    void print(etnte_node_ptr hdr, Context& c) {
      for (etnte_node_ptr etnte : HW::tree_inorder_range<etnte_algo>(hdr)) {
        directed_edge* de = Context::get_directed_edge(etnte);
        std::cout << fmt::format(" {}, {}", c.get_idx(de->v1), c.get_idx(Context::get_opposite(de)->v1));
      }
    }

    void print(et_node_ptr hdr, Context& c) {
      for (et_node_ptr et : HW::tree_inorder_range<et_algo>(hdr)) {
        if (!Context::is_loop_edge(et)) {
          directed_edge* de = Context::get_directed_edge(et);
          std::cout << fmt::format(" ({}, {})", c.get_idx(Context::get_opposite(de)->v1), c.get_idx(de->v1));  
        }
      }
    }

    // Returns true if reconnected.
    bool split_and_reconnect_tree_edge(directed_edge* ab) {
      auto ba = Context::get_opposite(ab);

      et_node_ptr et_ab = Context::get_tree_edge_entry(ab);
      et_node_ptr et_ba = Context::get_tree_edge_entry(ba);
      
      et_node_ptr et_hdr_a = et_algo::get_header(et_ab);
      etnte_node_ptr etnte_hdr_a = Context::get_etnte_header(et_hdr_a);

      et_node_ptr et_hdr_b = etPool_.acquire();
      etnte_node_ptr etnte_hdr_b = etntePool_.acquire();

      et_algo::init_header(et_hdr_b);
      etnte_algo::init_header(etnte_hdr_b);                           
      Context::set_etnte_header(et_hdr_b, etnte_hdr_b);                        

      etnte_context_operations<Context>::split(etnte_hdr_a, etnte_hdr_b, et_ab, et_ba);

      if (et_algo::less_than(et_ab, et_ba))
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ab, et_ba);
      else {
        cyclic<et_algo>::split(et_hdr_a, et_hdr_b, et_ba, et_ab);
        et_algo::swap_tree(et_hdr_a, et_hdr_b);        
      }      
      
                                                                
      directed_edge::set_null_et_entry(ab);
      directed_edge::set_null_et_entry(ba);      


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
      
      for (etnte_node_ptr etnte : HW::tree_inorder_range<etnte_algo>(etnte_hdr_a))
        if (etnte_algo::get_header(
          Context::get_non_tree_edge_entry(
            Context::get_opposite(
              Context::get_directed_edge(etnte)))) == etnte_hdr_b) {
          replacement_ab = etnte;
          break;
        }
     
      if (replacement_ab) {                                
        ab = Context::get_directed_edge(replacement_ab);
        ba = Context::get_opposite(ab);        
        etnte_node_ptr replacement_ba = Context::get_non_tree_edge_entry(ba);

        et_node_ptr et_rec_a = Context::get_ordering_vertex_entry(ab);
        et_node_ptr et_rec_b = Context::get_ordering_vertex_entry(ba);                                        

        Context::set_directed_edge(et_ab, ab);
        Context::set_directed_edge(et_ba, ba);

        Context::set_tree_edge_entry(ab, et_ab);
        Context::set_tree_edge_entry(ba, et_ba);

        // etnte
        cyclic<etnte_algo>::principal_cut(etnte_hdr_a, etnte_context_operations<Context>::lower_bound(etnte_hdr_a, et_rec_a));
        auto leastB = etnte_context_operations<Context>::lower_bound(etnte_hdr_b, et_rec_b);
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
      auto etnteHeader = Context::get_etnte_header(header);

      auto etnteNode = etnte_nt::get_left(etnteHeader);
      while (etnteNode != etnteHeader) {
        auto de = Context::get_directed_edge(etnteNode);
        directed_edge::set_null_et_entry(de);
        auto prev = etnteNode;
        etnteNode = etnte_algo::next_node(etnteNode);
        etntePool_.release(prev);
      }

      etntePool_.release(etnteHeader);


      auto etNode = et_nt::get_left(header);      
      while (etNode != header) {
        if (Context::is_loop_edge(etNode))
          Context::set_entry(Context::get_vertex(etNode), nullptr);
        else
          directed_edge::set_null_et_entry(Context::get_directed_edge(etNode));

        etNode = et_algo::next_node(etNode);
      }

      etNode = et_nt::get_left(header);
      while (etNode != header) {
        if (Context::is_loop_edge(etNode))
          visitor.vertex(Context::get_vertex(etNode));        

        auto prev = etNode;
        etNode = et_algo::next_node(etNode);
        etPool_.release(prev);
      }

      etPool_.release(header);
    }


    void remove_non_tree_edge(directed_edge* ab) {
      auto ba = Context::get_opposite(ab);
      auto etnteAB = Context::get_non_tree_edge_entry(ab);
      auto etnteBA = Context::get_non_tree_edge_entry(ba);
      auto etnteHeader = etnte_algo::get_header(etnteAB);
      
      etnte_algo::erase(etnteHeader, etnteAB);
      etnte_algo::erase(etnteHeader, etnteBA);

      directed_edge::set_null_et_entry(ab);
      directed_edge::set_null_et_entry(ba);

      etntePool_.release(etnteAB);
      etntePool_.release(etnteBA);
    }    
  };  
}
