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







    // void init_with_dc(dc_graph& g, dc_properties c) {
    //   using vertex_t = boost::graph_traits<dc_graph>::vertex_descriptor;
    //   using edge_t = boost::graph_traits<dc_graph>::edge_descriptor;
    //
    //   using dc_props = dc_properties;
    //
    //   std::cout << "starting build_euler_tour\n";
    //
    //   for (vertex_t v : g.range()) {
    //     et_ptr et_hdr = et_pool_.acquire();
    //     etnte_ptr etnte_hdr = etnte_pool_.acquire();
    //     et_algo::init_header(et_hdr);      
    //     etnte_algo::init_header(etnte_hdr);      
    //     dc_props::set_etnte_header(et_hdr, etnte_hdr);
    //
    //     et_ptr v_entry = et_pool_.acquire();
    //     dc_props::set_vertex(v_entry, v);
    //     dc_props::set_entry(v, v_entry);
    //     et_algo::push_back(et_hdr, v_entry); // TODO: can be optimised
    //   }
    //
    //   std::cout << "vertices initialised\n";
    //
    //   using namespace std::ranges;
    //   for (std::size_t i = 0; i < 2*g.edge_count(); ++i)
    //     if (edge_t ab = g.get_directed_edge(i);
    //       dc_props::is_null_entry(ab)) {
    //
    //     
    //
    //
    //       edge_t ba = ab->opposite;
    //
    //       et_ptr a_entry = dc_props::get_entry(ba->v1);
    //       et_ptr b_entry = dc_props::get_entry(ab->v1);
    //
    //       et_ptr a_hdr = et_algo::get_header(a_entry);
    //       et_ptr b_hdr = et_algo::get_header(b_entry);
    //
    //       if (a_hdr == b_hdr) { // connecting as non-tree edge
    //         // print(a_hdr, c);
    //         // print(b_hdr, c);
    //         // std::cout << '\n';
    //
    //               // etnte_ptr ab_entry = etnte_pool_.acquire();
    //               // dc_props::set_directed_edge(ab_entry, ab);
    //               // dc_props::set_non_tree_edge_entry(ab, ab_entry);
    //               //
    //               // etnte_ptr ba_entry = etnte_pool_.acquire();
    //               // dc_props::set_directed_edge(ba_entry, ba);
    //               // dc_props::set_non_tree_edge_entry(ba, ba_entry);
    //               //
    //               // if (etnte_ptr etnte_hdr = dc_props::get_etnte_header(a_hdr);
    //               //   etnte_algo::empty(etnte_hdr)) {
    //               //
    //               //   etnte_algo::push_back(etnte_hdr, ab_entry); // TODO: check order? probably not needed
    //               //   etnte_algo::push_back(etnte_hdr, ba_entry);
    //               // }
    //               // else {
    //               //   etnte_algo::insert_before(
    //               //     etnte_hdr,
    //               //     etnte_context_operations<dc_props>::lower_bound(etnte_hdr, a_entry),
    //               //     ab_entry);
    //               //
    //               //   etnte_algo::insert_before(
    //               //     etnte_hdr,
    //               //     etnte_context_operations<dc_props>::lower_bound(etnte_hdr, b_entry),
    //               //     ba_entry);
    //               // }
    //       }
    //       else { // connecting as tree edge
    //         // auto verify = [](etnte_ptr hdr) {
    //         //   
    //         //   if (auto* root = etnte_nt::get_parent(hdr); !(
    //         //     etnte_algo::verify(hdr) &&
    //         //     etnte_algo::calculate_subtree_size(root) == (root ? etnte_nt::get_size(root) : 0))) {
    //         //
    //         //     std::cout << "{A} INVALID\n";
    //         //     return false;
    //         //   }
    //         //
    //         //   return true;
    //         // };
    //
    //                     // { // etnte
    //                     //   etnte_ptr a_etnte_hdr = dc_props::get_etnte_header(a_hdr);
    //                     //   etnte_ptr b_etnte_hdr = dc_props::get_etnte_header(b_hdr);
    //                     //
    //                     //   // verify(a_etnte_hdr);
    //                     //   // verify(b_etnte_hdr);
    //                     //
    //                     //   if (!etnte_algo::empty(a_etnte_hdr)) {
    //                     //     etnte_ptr lb_a = etnte_context_operations<dc_props>::lower_bound(a_etnte_hdr, a_entry);
    //                     //     if (a_etnte_hdr != lb_a)
    //                     //       cyclic<etnte_algo>::cut(a_etnte_hdr, lb_a);
    //                     //     else {
    //                     //       // std::cout << "lb_a EQUAL HEADER 1\n";
    //                     //     }
    //                     //   }
    //                     //
    //                     //   etnte_ptr lb_b = nullptr;
    //                     //
    //                     //   if (!etnte_algo::empty(b_etnte_hdr)) {
    //                     //     lb_b = etnte_context_operations<dc_props>::lower_bound(b_etnte_hdr, b_entry);
    //                     //     if (b_etnte_hdr != lb_b) {
    //                     //       cyclic<etnte_algo>::cut_least_dropped(b_etnte_hdr, lb_b);
    //                     //       // cyclic<etnte_algo>::cut(b_etnte_hdr, lb_b);
    //                     //     }
    //                     //     else {
    //                     //       lb_b = etnte_nt::get_left(b_etnte_hdr);
    //                     //       etnte_algo::erase(b_etnte_hdr, lb_b);
    //                     //       // std::cout << "lb_b EQUAL HEADER 2\n";
    //                     //     }
    //                     //   }
    //                     //
    //                     //   if (etnte_algo::empty(b_etnte_hdr)) {
    //                     //     // nothing
    //                     //   }
    //                     //   else {
    //                     //     if (etnte_algo::empty(a_etnte_hdr)) {
    //                     //       etnte_algo::push_front(b_etnte_hdr, lb_b);
    //                     //       etnte_algo::swap_tree(a_etnte_hdr, b_etnte_hdr);
    //                     //     }
    //                     //     else {
    //                     //       // auto lb_b = etnte_nt::get_left(b_etnte_hdr);
    //                     //       // etnte_algo::erase(b_etnte_hdr, lb_b);
    //                     //       etnte_algo::join_trees(a_etnte_hdr, lb_b, b_etnte_hdr);
    //                     //     }
    //                     //   }
    //                     //
    //                     //
    //                     //   // if (!etnte_algo::empty(b_etnte_hdr)) {
    //                     //   //   if (etnte_algo::empty(a_etnte_hdr)) {
    //                     //   //
    //                     //   //     etnte_ptr lb_b = etnte_context_operations<dc_props>::lower_bound(b_etnte_hdr, b_entry);
    //                     //   //     if (b_etnte_hdr != lb_b) {
    //                     //   //       cyclic<etnte_algo>::cut(b_etnte_hdr, lb_b);
    //                     //   //     }
    //                     //   //     else {
    //                     //   //       std::cout << "lb_b EQUAL HEADER 0\n";
    //                     //   //     }
    //                     //   //
    //                     //   //     etnte_algo::swap_tree(a_etnte_hdr, b_etnte_hdr);
    //                     //   //
    //                     //   //
    //                     //   //   }
    //                     //   //   else {
    //                     //   //     etnte_ptr lb_a = etnte_context_operations<dc_props>::lower_bound(a_etnte_hdr, a_entry);
    //                     //   //     if (a_etnte_hdr != lb_a) {
    //                     //   //       cyclic<etnte_algo>::cut(a_etnte_hdr, lb_a);
    //                     //   //     }
    //                     //   //     else {
    //                     //   //       std::cout << "lb_a EQUAL HEADER 1\n";
    //                     //   //     }
    //                     //   //
    //                     //   //     etnte_ptr lb_b = etnte_context_operations<dc_props>::lower_bound(b_etnte_hdr, b_entry);
    //                     //   //     if (b_etnte_hdr != lb_b) {
    //                     //   //       cyclic<etnte_algo>::cut_least_dropped(a_etnte_hdr, lb_a);
    //                     //   //     }
    //                     //   //     else {
    //                     //   //       lb_b = etnte_nt::get_left(b_etnte_hdr);
    //                     //   //       etnte_algo::erase(b_etnte_hdr, lb_b);
    //                     //   //       std::cout << "lb_b EQUAL HEADER 2\n";
    //                     //   //     }
    //                     //   //
    //                     //   //     // cyclic<etnte_algo>::cut_least_dropped(b_etnte_hdr, lb_b);
    //                     //   //     etnte_algo::join_trees(a_etnte_hdr, lb_b, b_etnte_hdr);
    //                     //   //   }
    //                     //   // }
    //                     //   // else {
    //                     //   //   etnte_ptr lb_a = etnte_context_operations<dc_props>::lower_bound(a_etnte_hdr, a_entry);
    //                     //   //   if (a_etnte_hdr != lb_a) {
    //                     //   //     cyclic<etnte_algo>::cut(a_etnte_hdr, lb_a);
    //                     //   //   }
    //                     //   //   else {
    //                     //   //     std::cout << "lb_a EQUAL HEADER 4\n";
    //                     //   //   }
    //                     //   // }
    //                     //
    //                     //
    //                     //
    //                     //   // etnte_pool_.release(b_etnte_hdr);
    //                     // }
    //
    //         { // et
    //           et_ptr ab_entry = et_pool_.acquire();
    //           dc_props::set_directed_edge(ab_entry, ab);
    //           dc_props::set_tree_edge_entry(ab, ab_entry);
    //
    //           et_ptr ba_entry = et_pool_.acquire();
    //           dc_props::set_directed_edge(ba_entry, ba);
    //           dc_props::set_tree_edge_entry(ba, ba_entry);
    //
    //
    //           cyclic<et_algo>::cut(a_hdr, a_entry);
    //           cyclic<et_algo>::cut(b_hdr, b_entry);
    //
    //
    //
    //           et_algo::join_trees(a_hdr, ab_entry, b_hdr);
    //           et_algo::push_back(a_hdr, ba_entry);
    //
    //           // et_pool_.release(b_hdr);
    //         }
    //       }
    //     }
    //
    //   std::cout << "edges initialised\n";
    // }



    void init_with_dfs(dc_graph& g, dc_properties c) {           
      auto vertex_count = num_vertices(g);                 

      using vertex_t = dc_graph::vertex_descriptor;
      using edge_t = dc_graph::edge_descriptor;

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
        vertex_t v = g.get_vertex(i);

        if (color_map[v] != color_t::two_bit_black) {
          boost::depth_first_visit(g, v, visitor, color_map); // tree edges

          et_ptr et_hdr = et_algo::get_header(dc_properties::get_entry(v));
          etnte_ptr etnte_hdr = dc_properties::get_etnte_header(et_hdr);

          for (et_ptr et_entry : range<et_nt>(et_hdr))
            if (dc_properties::is_loop_edge(et_entry))
              for (edge_t de : g.range(g.get_idx(dc_properties::get_vertex(et_entry))))
                if (dc_properties::is_null_entry(de)) { // non-tree edges
                  etnte_ptr etnte_entry = etnte_pool_.acquire();
                  dc_properties::set_non_tree_edge_entry(de, etnte_entry); 
                  dc_properties::set_directed_edge(etnte_entry, de);

                  // avltree_algorithms_ext<etnte_traits>

                  // etnte_algo::insert_commit()



                  etnte_algo
                  ::push_back(etnte_hdr, etnte_entry);
                }
        }
      }
      // for (vertex* v : g.range())
      //   if (color_map[v] != color_t::two_bit_black) {
      //     boost::depth_first_visit(g, v, visitor, color_map);
      //   }
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

      et_ptr et_ab = Props::get_tree_edge_entry(ab);
      et_ptr et_ba = Props::get_tree_edge_entry(ba);
      
      et_ptr et_hdr_a = et_algo::get_header(et_ab);
      etnte_ptr etnte_hdr_a = Props::get_etnte_header(et_hdr_a);

      et_ptr et_hdr_b = et_pool_.acquire();
      etnte_ptr etnte_hdr_b = etnte_pool_.acquire();
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
        etnte_ptr root_a = etnte_nt::get_parent(etnte_hdr_a);
        etnte_ptr root_b = etnte_nt::get_parent(etnte_hdr_b);

        auto size_a = root_a ? etnte_nt::get_size(root_a) : 0;
        auto size_b = root_b ? etnte_nt::get_size(root_b) : 0;

        if (size_b < size_a) {
          std::swap(et_hdr_a, et_hdr_b);
          std::swap(etnte_hdr_a, etnte_hdr_b);
        }
      }

      etnte_ptr recon_ab_entry = nullptr;        
      
      for (etnte_ptr etnte : range<etnte_traits>(etnte_hdr_a))
        if (etnte_algo::get_header(
          Props::get_non_tree_edge_entry(
            Props::get_opposite(
              Props::get_directed_edge(etnte)))) == etnte_hdr_b) {
          recon_ab_entry = etnte;
          break;
        }
     
      if (recon_ab_entry) {                                
        ab = Props::get_directed_edge(recon_ab_entry);
        ba = Props::get_opposite(ab);        
        etnte_ptr replacement_ba_entry = Props::get_non_tree_edge_entry(ba);

        et_ptr recon_a_entry = Props::get_ordering_vertex_entry(ab);
        et_ptr recon_b_entry = Props::get_ordering_vertex_entry(ba);                                        

        Props::set_directed_edge(et_ab, ab);
        Props::set_directed_edge(et_ba, ba);

        Props::set_tree_edge_entry(ab, et_ab);
        Props::set_tree_edge_entry(ba, et_ba);

        // etnte
        cyclic<etnte_algo>::cut(etnte_hdr_a, etnte_context_operations<Props>::lower_bound(etnte_hdr_a, recon_a_entry));
        auto least_b = etnte_context_operations<Props>::lower_bound(etnte_hdr_b, recon_b_entry);
        cyclic<etnte_algo>::cut_least_dropped(etnte_hdr_b, least_b);
        etnte_algo::join_trees(etnte_hdr_a, least_b, etnte_hdr_b);

        etnte_algo::erase(etnte_hdr_a, recon_ab_entry);
        etnte_algo::erase(etnte_hdr_a, replacement_ba_entry);    

        // et
        cyclic<et_algo>::cut(et_hdr_a, recon_a_entry);
        cyclic<et_algo>::cut(et_hdr_b, recon_b_entry);                
                  
        et_algo::join_trees(et_hdr_a, et_ab, et_hdr_b);        
        et_algo::push_back(et_hdr_a, et_ba);                       
          
        return true;
      }      
     

      return false;
    }


        


    // template<class Visitor = component_visitor_empty>
    void remove_and_release_component(et_ptr header/*, Visitor& visitor = Visitor()*/) {            
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
        // if (Props::is_loop_edge(etNode))
        //   visitor.vertex(Props::get_vertex(etNode));        

        auto prev = etNode;
        etNode = et_algo::next_node(etNode);
        et_pool_.release(prev);
      }

      et_pool_.release(header);
    }


    
  };  
}
