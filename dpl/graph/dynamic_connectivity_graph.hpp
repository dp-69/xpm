#pragma once

#include "_OLD_general.hpp"
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/two_bit_color_map.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/trivial_value_traits.hpp>


#undef max

namespace HW::dynamic_connectivity
{
  // struct et_node;
  // struct etnte_node;
}

namespace HW::dynamic_connectivity
{
  struct dynamic_connectivity_graph;

  struct vertex;
  struct directed_edge;
}


namespace HW { namespace dynamic_connectivity
{
  // using et_node_ptr = HW::dynamic_connectivity::et_node*;
  // using etnte_node_ptr = etnte_node*;

  using vertex_ptr = vertex*;
  using directed_edge_ptr = directed_edge*;








  
  

  
  typedef default_list_node_traits<directed_edge> directed_edge_node_traits;



  typedef inorder_iter<directed_edge_node_traits> out_edge_iterator;
   

  struct directed_edge : non_copyable_movable
  {
    typedef directed_edge* node_ptr;    
//    typedef const directed_edge* const_node_ptr;

  private:
    friend struct default_list_node_traits<directed_edge>;

    // entry for the list of out edges from v0 pointing to v1
    node_ptr prev_;
    node_ptr next_;



//    typedef pointer_plus_aligned_least_significant_bits<void*, bool, 1> compression;
    typedef tagged_pointer_as_size_t<bool, 1> compression;

    size_t entry_type_; 

    

  public:
    vertex_ptr v1;  // pointing-in vertex  

    node_ptr opposite;

    
    static void set_null_et_entry(const node_ptr& x) {
      x->entry_type_ = 0;
    }

    static bool is_null_et_entry(const node_ptr& x) {
      return !x->entry_type_;
    }    



    static bool is_tree_edge(const node_ptr& x) {
      return compression::get_bits(x->entry_type_);
    }



    using et_node_ptr = et_traits::node_ptr;

    static et_node_ptr get_tree_edge_entry(const node_ptr& x) {
      return compression::get_pointer<et_node_ptr>(x->entry_type_);      
    }

    static void set_tree_edge_entry(const node_ptr& x, const et_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, true);      
    }



    static etnte_node_ptr get_non_tree_edge_entry(const node_ptr& x) {
      return compression::get_pointer<etnte_node_ptr>(x->entry_type_);      
    }

    static void set_non_tree_edge_entry(const node_ptr& x, const etnte_node_ptr& y) {
      compression::set_pointer_and_bits(x->entry_type_, y, false);      
    }    
  };

    
  typedef boost::intrusive::list<
    directed_edge,
    boost::intrusive::value_traits<boost::intrusive::trivial_value_traits<directed_edge_node_traits, boost::intrusive::normal_link>>,
    boost::intrusive::constant_time_size<false>> out_edge_list;
  


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
    friend pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_ptr, const dynamic_connectivity_graph&);

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

  typedef default_list_node_traits<vertex> vertex_node_traits;
 
  typedef boost::intrusive::list<
    vertex,
    boost::intrusive::value_traits<boost::intrusive::trivial_value_traits<vertex_node_traits, boost::intrusive::normal_link>>,
    boost::intrusive::constant_time_size<true>> vertex_list;

  typedef inorder_iter<vertex_node_traits> vertex_iterator;







    
  

  


    
  

    

  struct dynamic_connectivity_graph
  {
    typedef size_t vertices_size_type;
    typedef size_t edges_size_type;
    typedef size_t degree_size_type;

  private:
    vertex_list vertices_;

    friend void add_vertex(vertex_ptr v, dynamic_connectivity_graph& g);
    friend vertices_size_type num_vertices(const dynamic_connectivity_graph& g);
    friend pair<vertex_iterator, vertex_iterator> vertices(const dynamic_connectivity_graph& g);
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

  inline pair<vertex_iterator, vertex_iterator> vertices(const dynamic_connectivity_graph& g) {
	  return make_pair(vertex_iterator(g.vertices_.begin().pointed_node()), vertex_iterator(g.vertices_.end().pointed_node()));
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
  inline pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_ptr u, const dynamic_connectivity_graph& g) {
    return make_pair(out_edges_begin(u), out_edges_end(u));
  }
  
  // ReSharper disable once CppFunctionDoesntReturnValue
  inline vertex_ptr source(const directed_edge_ptr, const dynamic_connectivity_graph&) {}
  // ReSharper disable once CppFunctionDoesntReturnValue
  inline dynamic_connectivity_graph::degree_size_type out_degree(const vertex_ptr u, const dynamic_connectivity_graph&) {}
  //////////////////////////
  
 







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




  const vertex _NULL_VERTEX{};

  struct _generic_graph_traits {
  
  private:    
    struct traversal_tag :
      virtual boost::incidence_graph_tag
//      , virtual vertex_list_graph_tag
//      , virtual adjacency_graph_tag
//      , virtual edge_list_graph_tag
//      , virtual bidirectional_graph_tag 
    
    {};

  public:
    typedef vertex_ptr vertex_descriptor;
    typedef directed_edge_ptr edge_descriptor;
    
    typedef vertex_iterator vertex_iterator;
    typedef out_edge_iterator out_edge_iterator;
    
    
    typedef traversal_tag traversal_category;
    typedef boost::undirected_tag directed_category;

//    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;

    typedef dynamic_connectivity_graph::vertices_size_type vertices_size_type;
    typedef dynamic_connectivity_graph::vertices_size_type edges_size_type;
    typedef dynamic_connectivity_graph::degree_size_type degree_size_type;
        

    static vertex_descriptor null_vertex() {
      return const_cast<vertex_descriptor>(&_NULL_VERTEX);
    }
  };
}}



namespace boost {      
  template <class Graph>
  class empty_dfs_visitor : public default_dfs_visitor
  {      
    typedef Graph g;
    typedef typename graph_traits<Graph>::vertex_descriptor v;
    typedef typename graph_traits<Graph>::edge_descriptor e;

  public:      
    void initialize_vertex(v, const g&) {}
    void start_vertex(v, const g&) {}
    void discover_vertex(v, const g&) {}
    void finish_vertex(v, const g&) {}

    void examine_edge(e, const g&) {}
    void tree_edge(e, const g&) {}
    void back_edge(e, const g&) {}
    void forward_or_cross_edge(e, const g&) {}
  };         


  template <class Visitor1, class Visitor2>
  class merged_dfs_visitor 
  {      
    Visitor1 _v1;
    Visitor2 _v2;

  public:
    merged_dfs_visitor(Visitor1& v1, Visitor2& v2)
      : _v1(v1), _v2(v2) {}

    template<class Vertex, class Graph>
    void initialize_vertex(Vertex v, const Graph& g) {
      _v1.initialize_vertex(v, g);
      _v2.initialize_vertex(v, g);
    }

    template<class Vertex, class Graph>
    void start_vertex(Vertex v, const Graph& g) {
      _v1.start_vertex(v, g);
      _v2.start_vertex(v, g);  
    }

    template<class Vertex, class Graph>
    void discover_vertex(Vertex v, const Graph& g) {
      _v1.discover_vertex(v, g);
      _v2.discover_vertex(v, g);
    }

    template<class Vertex, class Graph>
    void finish_vertex(Vertex v, const Graph& g) {
      _v1.finish_vertex(v, g);
      _v2.finish_vertex(v, g);
    }

    template<class Edge, class Graph>
    void examine_edge(Edge e, const Graph& g) {
      _v1.examine_edge(e, g);
      _v2.examine_edge(e, g);
    }

    template<class Edge, class Graph>
    void tree_edge(Edge e, const Graph& g) {
      _v1.tree_edge(e, g);
      _v2.tree_edge(e, g);
    }

    template<class Edge, class Graph>
    void back_edge(Edge e, const Graph& g) {
      _v1.back_edge(e, g);
      _v2.back_edge(e, g);
    }

    template<class Edge, class Graph>
    void forward_or_cross_edge(Edge e, const Graph& g) {
      _v1.forward_or_cross_edge(e, g);
      _v2.forward_or_cross_edge(e, g);
    }
  };         

  template<class Visitor1, class Visitor2>
  merged_dfs_visitor<Visitor1, Visitor2> merge_visitors(Visitor1 v1, Visitor2 v2) {
    return merged_dfs_visitor<Visitor1, Visitor2>(v1, v2);
  }


  template<>
  struct graph_traits<HW::dynamic_connectivity::dynamic_connectivity_graph> :
           HW::dynamic_connectivity::_generic_graph_traits {};
}