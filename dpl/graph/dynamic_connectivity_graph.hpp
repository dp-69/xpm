#pragma once

#include "_OLD_general.hpp"
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/two_bit_color_map.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/trivial_value_traits.hpp>


#undef max


namespace HW::dynamic_connectivity
{
  struct vertex;

  using vertex_node_traits = default_list_node_traits<vertex>;
  using vertex_iterator = inorder_iter<vertex_node_traits>;
  using vertex_list = boost::intrusive::list<
    vertex,
    boost::intrusive::value_traits<boost::intrusive::trivial_value_traits<vertex_node_traits>>,
    boost::intrusive::constant_time_size<true>>;


  struct directed_edge;

  using directed_edge_node_traits = default_list_node_traits<directed_edge>;
  using out_edge_iterator = inorder_iter<directed_edge_node_traits>;
  using out_edge_list = boost::intrusive::list<
    directed_edge,
    boost::intrusive::value_traits<boost::intrusive::trivial_value_traits<directed_edge_node_traits>>,
    boost::intrusive::constant_time_size<false>>;

  struct dynamic_connectivity_graph;
}







namespace HW::dynamic_connectivity
{
  using vertex_ptr = vertex*;
  using directed_edge_ptr = directed_edge*;

  
 
   

  struct directed_edge : non_copyable_movable
  {
  private:
    friend struct default_list_node_traits<directed_edge>;

    // entry for the list of out edges from v0 pointing to v1
    directed_edge_ptr prev_;
    directed_edge_ptr next_;



//    typedef pointer_plus_aligned_least_significant_bits<void*, bool, 1> compression;
    typedef tagged_pointer_as_size_t<bool, 1> compression;

    size_t entry_type_; 

    

  public:
    vertex_ptr v1;  // pointing-in vertex  

    directed_edge_ptr opposite;

    
    static void set_null_et_entry(const directed_edge_ptr& x) {
      x->entry_type_ = 0;
    }

    static bool is_null_et_entry(const directed_edge_ptr& x) {
      return !x->entry_type_;
    }    









    static bool is_tree_edge(const directed_edge_ptr& x) {
      return compression::get_bits(x->entry_type_);
    }






    using et_node_ptr = et_traits::node_ptr;

    static et_node_ptr get_tree_edge_entry(const directed_edge_ptr& x) {
      return compression::get_pointer<et_node_ptr>(x->entry_type_);      
    }

    static void set_tree_edge_entry(const directed_edge_ptr& x, const et_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, true);      
    }

    using etnte_node_ptr = etnte_traits::node_ptr;

    static etnte_node_ptr get_non_tree_edge_entry(const directed_edge_ptr& x) {
      return compression::get_pointer<etnte_node_ptr>(x->entry_type_);      
    }

    static void set_non_tree_edge_entry(const directed_edge_ptr& x, const etnte_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, false);      
    }    
  };

    
  
  


  struct vertex : non_copyable_movable
  {
    typedef vertex* node_ptr;

  private:    
    friend struct default_list_node_traits<vertex>;

    // intrusive list of graph vertices
    node_ptr prev_;
    node_ptr next_;          
    
    out_edge_list out_edges_;

    
    friend void clear_out_edges(vertex_ptr, dynamic_connectivity_graph&);
    friend void add_edge(vertex&, vertex&, directed_edge&, directed_edge&, dynamic_connectivity_graph&);
    friend std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_ptr, const dynamic_connectivity_graph&);

    friend out_edge_iterator out_edges_begin(const vertex_ptr);
    friend out_edge_iterator out_edges_end(const vertex_ptr);

//    friend class euler_tour_dynamic_connectivity_context;
//    friend struct euler_tour_non_tree_edge_node_traits;
//    friend class euler_tour_visitor;
//    friend class pnm_et_dc_context;
    

//    friend out_edge_iterator out_edges_begin(const vertex_ptr, const dc_graph&);
//    friend out_edge_iterator out_edges_end(const vertex_ptr, const dc_graph&);
  public:    
    size_t row_idx_ = 0;  // relative index, not unique identifier    
//    static const size_t invalid_rel_idx = ~size_t(0) >> 2;


    et_traits::node_ptr et_entry_;
    bool visited = false;
//    bool trapped = false;
  };

  







    
  

  
  






































































    
  

    

  struct dynamic_connectivity_graph
  {
    typedef size_t vertices_size_type;
    typedef size_t edges_size_type;
    typedef size_t degree_size_type;

  private:
    vertex_list vertices_;

    friend void add_vertex(vertex_ptr v, dynamic_connectivity_graph& g);
    friend vertices_size_type num_vertices(const dynamic_connectivity_graph& g);
    friend std::pair<vertex_iterator, vertex_iterator> vertices(const dynamic_connectivity_graph& g);
  };

  inline void add_vertex(vertex_ptr v, dynamic_connectivity_graph& g) {            
    g.vertices_.push_back(*v);    
  }

  inline void remove_vertex(vertex_ptr v, dynamic_connectivity_graph&) {
    vertex_list::node_algorithms::unlink(v);            
  }

  inline void clear_out_edges(vertex_ptr v, dynamic_connectivity_graph&) {
    v->out_edges_.clear();
  }

  inline void add_edge(vertex& v, vertex& u, directed_edge& vu, directed_edge& uv, dynamic_connectivity_graph&) {            
    vu.v1 = &u;
    v.out_edges_.push_back(vu);
    
    uv.v1 = &v;
    u.out_edges_.push_back(uv);

    vu.opposite = &uv;
    uv.opposite = &vu;
  }

  inline void remove_edge(const directed_edge_ptr& vu, dynamic_connectivity_graph&) {        
    out_edge_list::node_algorithms::unlink(vu->opposite);
    out_edge_list::node_algorithms::unlink(vu);
  }

  inline dynamic_connectivity_graph::vertices_size_type num_vertices(const dynamic_connectivity_graph& g) {    
    return g.vertices_.size();
  }

  inline std::pair<vertex_iterator, vertex_iterator> vertices(const dynamic_connectivity_graph& g) {
	  return std::make_pair(vertex_iterator(g.vertices_.begin().pointed_node()), vertex_iterator(g.vertices_.end().pointed_node()));
  }    
   
  inline out_edge_iterator out_edges_begin(const vertex_ptr v) {
    return out_edge_iterator(v->out_edges_.begin().pointed_node());    
  }

  inline out_edge_iterator out_edges_end(const vertex_ptr v) {
    return out_edge_iterator(v->out_edges_.end().pointed_node());    
  }

  struct out_edges_range
  {
    const vertex_ptr v;
    explicit out_edges_range(const vertex_ptr v) : v(v) {}
    out_edge_iterator begin() const { return out_edges_begin(v); }
    out_edge_iterator end() const { return out_edges_end(v); }
  };

  
  inline vertex_ptr target(const directed_edge_ptr e, const dynamic_connectivity_graph&) {
    return e->v1;
  }



  // for Boost concept check
  inline std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_ptr u, const dynamic_connectivity_graph& g) {
    return std::make_pair(out_edges_begin(u), out_edges_end(u));
  }
  
  inline vertex_ptr source(const directed_edge_ptr, const dynamic_connectivity_graph&) {}
  inline dynamic_connectivity_graph::degree_size_type out_degree(const vertex_ptr u, const dynamic_connectivity_graph&) {}
  
 







  struct vertex_color_property_map
  {       
    typedef boost::read_write_property_map_tag category;
    typedef boost::two_bit_color_type value_type; 
    typedef value_type& reference;
    typedef vertex_ptr key_type;
    

//  private:
//    typedef pointer_plus_aligned_least_significant_bits<euler_tour_node_ptr, value_type, 2> compression;
    typedef intergral_plus_shifted_least_significant_bits<size_t, value_type, 2> compression;

    static size_t& field(const vertex_ptr v) { return v->row_idx_; }

  public:                
    static void compress_and_init(const vertex_ptr v) {      
      set_value(v, field(v));
      set_color(v, boost::color_traits<value_type>::white());
    }

    static void decompress_and_finish(const vertex_ptr v) {
      field(v) = compression::get_value(field(v));
    }

    static value_type get_color(const vertex_ptr v) {
      return compression::get_bits(field(v));   
    }

    static void set_color(const vertex_ptr v, value_type color) {
      return compression::set_bits(field(v), color);   
    }

    static size_t get_value(const vertex_ptr v) {
      return compression::get_value(field(v));
    }

    static void set_value(const vertex_ptr v, size_t value) {
      compression::set_value(field(v), value);
    }
  };

  inline vertex_color_property_map::value_type get(const vertex_color_property_map&, const vertex_ptr v) {
    return vertex_color_property_map::get_color(v);    
  }

  inline void put(const vertex_color_property_map&, const vertex_ptr v, vertex_color_property_map::value_type color) {
    vertex_color_property_map::set_color(v, color);      
  }



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


namespace boost
{
  template <>
  struct graph_traits<HW::dynamic_connectivity::dynamic_connectivity_graph>
  {
  private:
    struct traversal_tag : virtual boost::incidence_graph_tag {};

  public:
    using vertex_descriptor = HW::dynamic_connectivity::vertex_ptr;
    using edge_descriptor = HW::dynamic_connectivity::directed_edge_ptr;

    using vertex_iterator = HW::dynamic_connectivity::vertex_iterator;
    using out_edge_iterator = HW::dynamic_connectivity::out_edge_iterator;


    using traversal_category = traversal_tag;
    using directed_category = boost::undirected_tag;

    //    typedef allow_parallel_edge_tag edge_parallel_category;
    using edge_parallel_category = boost::disallow_parallel_edge_tag;

    using vertices_size_type = HW::dynamic_connectivity::dynamic_connectivity_graph::vertices_size_type;
    using edges_size_type = HW::dynamic_connectivity::dynamic_connectivity_graph::vertices_size_type;
    using degree_size_type = HW::dynamic_connectivity::dynamic_connectivity_graph::degree_size_type;

    static vertex_descriptor null_vertex() {
      return const_cast<vertex_descriptor>(&HW::dynamic_connectivity::_NULL_VERTEX);
    }
  };
}
