﻿#pragma once

#include "_OLD_general.hpp"
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/trivial_value_traits.hpp>

namespace dpl::graph
{
  struct vertex;

  using vertex_const_iterator = std::vector<vertex>::const_iterator;
  using vertex_iterator = std::vector<vertex>::iterator;

  struct directed_edge;

  using directed_edge_node_traits = HW::default_list_node_traits<directed_edge>;
  using out_edge_iterator = HW::inorder_iter<directed_edge_node_traits>;
  using out_edge_list = boost::intrusive::list<
    directed_edge,
    boost::intrusive::value_traits<boost::intrusive::trivial_value_traits<directed_edge_node_traits>>,
    boost::intrusive::constant_time_size<false>>;

  class dc_graph;
}



namespace dpl::graph
{
  struct directed_edge
  {
  private:
    friend struct HW::default_list_node_traits<directed_edge>;
    friend class etnte_context;
    friend void add_edge(vertex&, vertex&, directed_edge&, directed_edge&, dc_graph&);
    friend void remove_edge(directed_edge* vu, dc_graph&);

    directed_edge* prev_;
    directed_edge* next_;

    directed_edge* opposite_;

    size_t entry_type_; 
  public:

    vertex* v1;  // pointing-to vertex  

    static void set_null_et_entry(directed_edge* x) {
      x->entry_type_ = 0;
    }

    static bool is_null_et_entry(const directed_edge* x) {
      return !x->entry_type_;
    }    
  };

  struct vertex
  {
  private:    
    friend class etnte_context;
    
    out_edge_list out_edges_;
    
    friend void clear_out_edges(vertex*, dc_graph&);
    friend void add_edge(vertex&, vertex&, directed_edge&, directed_edge&, dc_graph&);
    friend std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex*, const dc_graph&);

    friend out_edge_iterator out_edges_begin(const vertex*);
    friend out_edge_iterator out_edges_end(const vertex*);

    et_traits::node_ptr et_entry_;
  };

    

  class dc_graph
  {
    std::vector<vertex> vertices_;
    // friend void add_vertex(vertex* v, dc_graph&);

  public:



    using vertices_size_type = size_t;
    using edges_size_type = size_t;
    using degree_size_type = size_t;


    friend vertices_size_type num_vertices(const dc_graph&);
    friend vertex_const_iterator begin(const dc_graph&);
    friend vertex_const_iterator end(const dc_graph&);

    friend vertex_iterator begin(dc_graph&);
    friend vertex_iterator end(dc_graph&);

    auto& vertices() { return vertices_; }
  };

  // inline void add_vertex(vertex* v, dc_graph& g) {            
  //   g.vertices_.push_back(*v);    
  // }

  // inline void remove_vertex(vertex* v, dc_graph&) {
  //   vertex_list::node_algorithms::unlink(v);            
  // }

  inline void clear_out_edges(vertex* v, dc_graph&) {
    v->out_edges_.clear();
  }

  inline void add_edge(vertex& v, vertex& u, directed_edge& vu, directed_edge& uv, dc_graph&) {            
    vu.v1 = &u;
    v.out_edges_.push_back(vu);
    
    uv.v1 = &v;
    u.out_edges_.push_back(uv);

    vu.opposite_ = &uv;
    uv.opposite_ = &vu;
  }

  inline void remove_edge(directed_edge* vu, dc_graph&) {
    out_edge_list::node_algorithms::unlink(vu->opposite_);
    out_edge_list::node_algorithms::unlink(vu);
  }

  inline dc_graph::vertices_size_type num_vertices(const dc_graph& g) {    
    return g.vertices_.size();
  }

  inline vertex_const_iterator begin(const dc_graph& g) { return g.vertices_.begin(); }
  inline vertex_const_iterator end(const dc_graph& g) { return g.vertices_.end(); }

  inline vertex_iterator begin(dc_graph& g) { return g.vertices_.begin(); }
  inline vertex_iterator end(dc_graph& g) { return g.vertices_.end(); }

  inline std::pair<vertex_const_iterator, vertex_const_iterator> vertices(const dc_graph& g) {
	  return std::make_pair(begin(g), end(g));
  }

  inline auto range(const dc_graph& g) { return std::ranges::subrange{begin(g), end(g)}; }
  inline auto range(dc_graph& g) { return std::ranges::subrange{begin(g), end(g)}; }    
   
  inline out_edge_iterator out_edges_begin(const vertex* v) {
    return out_edge_iterator(v->out_edges_.begin().pointed_node());    
  }

  inline out_edge_iterator out_edges_end(const vertex* v) {
    return out_edge_iterator(v->out_edges_.end().pointed_node());    
  }

  // struct out_edges_range
  // {
  //   const vertex_ptr v;
  //   explicit out_edges_range(const vertex_ptr v) : v(v) {}
  //   out_edge_iterator begin() const { return out_edges_begin(v); }
  //   out_edge_iterator end() const { return out_edges_end(v); }
  // };

  
  inline vertex* target(const directed_edge* e, const dc_graph&) {
    return e->v1;
  }



  // for Boost concept check
  inline std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex* u, const dc_graph& g) {
    return std::make_pair(out_edges_begin(u), out_edges_end(u));
  }
  
  inline vertex* source(const directed_edge*, const dc_graph&) {}
  inline dc_graph::degree_size_type out_degree(const vertex* u, const dc_graph&) {}
  
 







  // struct vertex_color_property_map
  // {
  //   using category = boost::read_write_property_map_tag;
  //   using value_type = boost::two_bit_color_type;
  //   using reference = value_type&;
  //   using key_type = vertex_ptr;
  //   
  //   using compression = HW::intergral_plus_shifted_least_significant_bits<size_t, value_type, 2>;
  //
  //   static size_t& field(const vertex_ptr v) { return v->row_idx_; }
  //
  //   static void compress_and_init(const vertex_ptr v) {
  //     compression::set_value(field(v), field(v));
  //     set_color(v, boost::color_traits<value_type>::white());
  //   }
  //
  //   static void decompress_and_finish(const vertex_ptr v) {
  //     field(v) = compression::get_value(field(v));
  //   }
  //
  //   static value_type get_color(const vertex_ptr v) {
  //     // return v->color;
  //     return compression::get_bits(field(v));   
  //   }
  //
  //   static void set_color(const vertex_ptr v, value_type color) {
  //     // v->color = color;
  //     return compression::set_bits(field(v), color);   
  //   }
  // };
  //
  // inline vertex_color_property_map::value_type get(const vertex_color_property_map&, const vertex_ptr v) {
  //   return vertex_color_property_map::get_color(v);    
  // }
  //
  // inline void put(const vertex_color_property_map&, const vertex_ptr v, vertex_color_property_map::value_type color) {
  //   vertex_color_property_map::set_color(v, color);      
  // }



//  struct vertex_component_property_map : b::put_get_helper<decremental_connectivity_graph::vertices_size_type&, vertex_component_property_map>
//  {    
//    typedef b::writable_property_map_tag category;
//    typedef decremental_connectivity_graph::vertices_size_type value_type; 
//    typedef value_type& reference;
//    typedef vertex_ptr key_type;
//
//    value_type* componentIdx_;
//
//    reference operator[](key_type v) const {                        
//      return componentIdx_[vertex_color_property_map::compression::get_value(v->rel_idx_)];
//    }
//  };



//  inline vertex_component_property_map::value_type get(const vertex_component_property_map& pm, const graph::vertex_descriptor v) {
//    return pm[v];    
//  }
//
//  inline void put(const vertex_component_property_map& pm, const graph::vertex_descriptor v, vertex_component_property_map::value_type componentIdx) {
//    pm[v] = componentIdx;
//  }    




  static const vertex _NULL_VERTEX = {};

}


template <>
struct boost::graph_traits<dpl::graph::dc_graph>
{
private:
  struct traversal_tag : virtual incidence_graph_tag {};

public:
  using vertex_descriptor = dpl::graph::vertex*;
  using edge_descriptor = dpl::graph::directed_edge*;

  using out_edge_iterator = dpl::graph::out_edge_iterator;

  using traversal_category = traversal_tag;
  using directed_category = undirected_tag;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = dpl::graph::dc_graph::vertices_size_type;
  using edges_size_type = dpl::graph::dc_graph::vertices_size_type;
  using degree_size_type = dpl::graph::dc_graph::degree_size_type;

  static vertex_descriptor null_vertex() {
    return const_cast<vertex_descriptor>(&dpl::graph::_NULL_VERTEX);
  }
};
