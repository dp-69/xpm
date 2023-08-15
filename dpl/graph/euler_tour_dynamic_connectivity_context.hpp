#pragma once

#include "euler_tour_visitor.hpp"
// #include "HW/pore_network_modelling/row_idx_populate_visitor.h" //TODO

namespace HW { namespace dynamic_connectivity
{      
  using namespace boost;

  struct component_visitor_empty
  {
    void vertex(vertex_ptr) {}
    void directed_edge(directed_edge_ptr) {}
  };
    
  struct component_visitor_count
  {
    size_t _vertices = 0;
    size_t _directedEdges = 0;

    void vertex(vertex_ptr) { ++_vertices; }
    void directed_edge(directed_edge_ptr) { ++_directedEdges; }

    size_t vertices() { return _vertices; }
    size_t edges() { return _directedEdges/2; }
  };


  struct dc_stat
  {
    size_t nteCountOther;
    size_t nteCountActive;
    size_t nteReconnectingFirst;
  };

  class euler_tour_dynamic_connectivity_context : non_copyable_movable
  {
    using et_node_ptr = et_traits::node_ptr;
    using etnte_node_ptr = etnte_traits::node_ptr;

    const static auto stack_capacity = 256;

    using et_nt = et_traits;
    using etnte_nt = etnte_traits;
    
  public:    
  
  #ifdef CONNECTIVITY_DIAGNOSTICS
  vector<dc_stat> dcStatistics;
  #endif 

    smart_pool<et_traits::node> etPool_;
    smart_pool<etnte_traits::node> etntePool_;


    std::vector<et_node_ptr> initial_components_;
   
    template<class Visitor>
    void init(const dynamic_connectivity_graph& g, vertex_ptr rootVertex, Visitor v) {            
      auto vertexCount = num_vertices(g);                 
      auto treeEdgeStack = std::vector<et_node_ptr>(vertexCount + 1);     
          
      execute_dfs(g, rootVertex,              
        merge_visitors(euler_tour_visitor(etPool_, etntePool_, treeEdgeStack.data(), initial_components_), v));
    }
    
    void init(const dynamic_connectivity_graph& g, vertex_ptr rootVertex) {           
      init(g, rootVertex, empty_dfs_visitor<dynamic_connectivity_graph>());
    }

    // Returns true if reconnected.
    bool split_and_reconnect_tree_edge(directed_edge_ptr ab) {
      auto ba = ab->opposite;

      auto etAB = directed_edge::get_tree_edge_entry(ab);
      auto etBA = directed_edge::get_tree_edge_entry(ba);
      
      auto etHeaderA = et_algo::get_header(etAB);
      auto etnteHeaderA = etnte_context::get_non_tree_edge_header(etHeaderA);

      auto etHeaderB = etPool_.acquire();
      et_algo::init_header(etHeaderB);
        
      auto etnteHeaderB = etntePool_.acquire();
      etnte_algo::init_header(etnteHeaderB);                           
      etnte_context::set_non_tree_edge_header(etHeaderB, etnteHeaderB);                        

      etnte_context_operations<etnte_context>::split(etnteHeaderA, etnteHeaderB, etAB, etBA);

      if (et_algo::less_than(etAB, etBA))
        dpl::graph::cyclic<et_algo>::split(etHeaderA, etHeaderB, etAB, etBA);
      else {
        dpl::graph::cyclic<et_algo>::split(etHeaderA, etHeaderB, etBA, etAB);
        et_algo::swap_tree(etHeaderA, etHeaderB);        
      }      
      
                                                                
      directed_edge::set_null_et_entry(ab);
      directed_edge::set_null_et_entry(ba);      
              
      auto etnteSizeA = etnte_nt::get_size(etnte_nt::get_parent(etnteHeaderA));
      auto etnteSizeB = etnte_nt::get_size(etnte_nt::get_parent(etnteHeaderB));

      if (etnteSizeB < etnteSizeA) {
        std::swap(etHeaderA, etHeaderB);
        std::swap(etnteHeaderA, etnteHeaderB);
      }

      etnte_node_ptr replacementAB = nullptr;        

      #ifdef CONNECTIVITY_DIAGNOSTICS

      size_t i = 0;                                

      for (const auto& etnte : tree_inorder_range<etnte_algo>(etnteHeaderA)) {
        ++i;

        if (etnte_algo::get_header(directed_edge::get_non_tree_edge_entry(etnte_nt::get_directed_edge(etnte)->opposite)) == etnteHeaderB) {            
          dcStatistics.push_back(dc_stat{  
            etnte_nt::get_size(etnte_nt::get_parent(etnteHeaderB)),
            etnte_nt::get_size(etnte_nt::get_parent(etnteHeaderA)),            
//            et_nt::get_size(et_nt::get_parent(etHeaderB)),
//            et_nt::get_size(et_nt::get_parent(etHeaderA)),      
            i});

          replacementAB = etnte;
          break;
        }
      }
       
      #else
      for (const auto& etnte : tree_inorder_range<etnte_algo>(etnteHeaderA))
        if (etnte_algo::get_header(directed_edge::get_non_tree_edge_entry(etnte_context::get_directed_edge(etnte)->opposite)) == etnteHeaderB) {
          replacementAB = etnte;
          break;
        }
      #endif

     
      if (replacementAB) {                                
        ab = etnte_context::get_directed_edge(replacementAB);
        ba = ab->opposite;        
        auto replacementBA = directed_edge::get_non_tree_edge_entry(ba);

        auto etRecA = etnte_context::get_vertex_entry(ab);
        auto etRecB = etnte_context::get_vertex_entry(ba);                                        

        etnte_context::set_directed_edge(etAB, ab);
        etnte_context::set_directed_edge(etBA, ba);

        directed_edge::set_tree_edge_entry(ab, etAB);
        directed_edge::set_tree_edge_entry(ba, etBA);
                

        // etnte
        dpl::graph::cyclic<etnte_algo>::principal_cut(etnteHeaderA, etnte_context_operations<etnte_context>::lower_bound(etnteHeaderA, etRecA));
        auto leastB = etnte_context_operations<etnte_context>::lower_bound(etnteHeaderB, etRecB);
        dpl::graph::cyclic<etnte_algo>::principal_cut_least_dropped(etnteHeaderB, leastB);
        etnte_algo::join_trees(etnteHeaderA, leastB, etnteHeaderB);

        etnte_algo::erase(etnteHeaderA, replacementAB);
        etnte_algo::erase(etnteHeaderA, replacementBA);    


        // et
        dpl::graph::cyclic<et_algo>::principal_cut(etHeaderA, etRecA);
        dpl::graph::cyclic<et_algo>::principal_cut(etHeaderB, etRecB);                
                  
        et_algo::join_trees(etHeaderA, etAB, etHeaderB);        
        et_algo::push_back(etHeaderA, etBA);                       
          
        return true;
      }      
     

      return false;
    }


        


    template<class Visitor = component_visitor_empty>
    void remove_and_release_component(et_node_ptr header, Visitor& visitor = Visitor()) {            
      auto etnteHeader = etnte_context::get_non_tree_edge_header(header);

      auto etnteNode = etnte_nt::get_left(etnteHeader);
      while (etnteNode != etnteHeader) {
        auto de = etnte_context::get_directed_edge(etnteNode);
        directed_edge::set_null_et_entry(de);
        auto prev = etnteNode;
        etnteNode = etnte_algo::next_node(etnteNode);
        etntePool_.release(prev);
      }

      etntePool_.release(etnteHeader);


      auto etNode = et_nt::get_left(header);      
      while (etNode != header) {
        if (etnte_context::is_loop_edge(etNode))
          etnte_context::get_vertex(etNode)->et_entry_ = nullptr;
        else
          directed_edge::set_null_et_entry(etnte_context::get_directed_edge(etNode));

        etNode = et_algo::next_node(etNode);
      }

      etNode = et_nt::get_left(header);
      while (etNode != header) {
        if (etnte_context::is_loop_edge(etNode))
          visitor.vertex(etnte_context::get_vertex(etNode));        

        auto prev = etNode;
        etNode = et_algo::next_node(etNode);
        etPool_.release(prev);
      }

      etPool_.release(header);
    }


    void remove_non_tree_edge(directed_edge_ptr ab) {
      auto ba = ab->opposite;
      auto etnteAB = directed_edge::get_non_tree_edge_entry(ab);
      auto etnteBA = directed_edge::get_non_tree_edge_entry(ba);
      auto etnteHeader = etnte_algo::get_header(etnteAB);
      
      etnte_algo::erase(etnteHeader, etnteAB);
      etnte_algo::erase(etnteHeader, etnteBA);

      directed_edge::set_null_et_entry(ab);
      directed_edge::set_null_et_entry(ba);

      etntePool_.release(etnteAB);
      etntePool_.release(etnteBA);
    }    
  };  
}}
