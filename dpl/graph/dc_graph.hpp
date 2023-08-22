#pragma once

namespace dpl::graph
{
  template<typename Pointer>
  struct iterator_deference_pointer
  {
    Pointer ptr;

    using iterator_category = std::forward_iterator_tag;
    using value_type = Pointer;
    using difference_type = std::ptrdiff_t;
    using distance_type = std::ptrdiff_t;
    using pointer = Pointer;
    using reference = Pointer;

    pointer operator*() const { return ptr; }
    pointer operator->() { return ptr; }

    iterator_deference_pointer& operator++() {
      ++ptr;
      return *this;
    }  

    iterator_deference_pointer operator++(int) {
      iterator_deference_pointer temp = *this;
      ++*this;
      return temp;
    }

    friend bool operator== (const iterator_deference_pointer& a, const iterator_deference_pointer& b) { return a.ptr == b.ptr; }
    friend bool operator!= (const iterator_deference_pointer& a, const iterator_deference_pointer& b) { return a.ptr != b.ptr; }
  };
}



namespace dpl::graph
{
  class vertex
  {
    friend class etnte_properties;

    et_traits::node_ptr entry_;

  public:

    vertex() = default;
    vertex(const vertex& other) = delete;
    vertex(vertex&& other) noexcept = delete;
    vertex& operator=(const vertex& other) = delete;
    vertex& operator=(vertex&& other) noexcept = delete;
  };

  class directed_edge
  {
    friend class etnte_properties;

    size_t entry_;

  public:
    directed_edge() = default;
    directed_edge(const directed_edge& other) = delete;
    directed_edge(directed_edge&& other) noexcept = delete;
    directed_edge& operator=(const directed_edge& other) = delete;
    directed_edge& operator=(directed_edge&& other) noexcept = delete;

    directed_edge* opposite;
    vertex* v1;  // pointing-to vertex  
  };


  using vertex_iter_custom = iterator_deference_pointer<vertex*>;
  using out_edge_iterator = iterator_deference_pointer<directed_edge*>;

  class dc_graph
  {
    std::vector<vertex> vertices_;
    std::vector<directed_edge> directed_edges_;
    std::unique_ptr<std::size_t[]> directed_edges_end_;

    // friend std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex*, const dc_graph&);

  public:
    using vertices_size_type = std::size_t;
    using edges_size_type = std::size_t;
    using degree_size_type = std::size_t;

    explicit dc_graph(vertices_size_type count)
      : vertices_(count) {}

    dc_graph(const dc_graph& other) = delete;
    dc_graph(dc_graph&& other) noexcept = delete;
    dc_graph& operator=(const dc_graph& other) = delete;
    dc_graph& operator=(dc_graph&& other) noexcept = delete;

    const auto& vertices() const { return vertices_; }
    auto& vertices() { return vertices_; }

    const auto& directed_edges() const { return directed_edges_; }
    auto& directed_edges() { return directed_edges_; }

    const auto& directed_edges_end() const { return directed_edges_end_; }
    auto& directed_edges_end() { return directed_edges_end_; }
  };

  inline dc_graph::vertices_size_type num_vertices(const dc_graph& g) {    
    return g.vertices().size();
  }

  inline vertex_iter_custom begin(dc_graph& g) { return {g.vertices().data()}; }
  inline vertex_iter_custom end(dc_graph& g) { return {g.vertices().data() + g.vertices().size()}; }

  inline std::pair<vertex_iter_custom, vertex_iter_custom> vertices(dc_graph& g) {
	  return std::make_pair(begin(g), end(g));
  }

  inline vertex* target(const directed_edge* e, const dc_graph&) {
    return e->v1;
  }

  inline std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex* u, const dc_graph& g) {
    auto idx = u - g.vertices().data();

    return std::make_pair(
      out_edge_iterator{
        const_cast<directed_edge*>(g.directed_edges().data() + (idx == 0 ? 0 : g.directed_edges_end()[idx - 1]))},
      out_edge_iterator{
        const_cast<directed_edge*>(g.directed_edges().data() + g.directed_edges_end()[idx])}
      );
  }
  
  inline vertex* source(const directed_edge*, const dc_graph&) {}
  inline dc_graph::degree_size_type out_degree(const vertex* u, const dc_graph&) {}


  inline auto range(dc_graph& g) {
    return std::ranges::subrange{begin(g), end(g)};
  }    
}


template <>
struct boost::graph_traits<dpl::graph::dc_graph>
{
  using vertex_descriptor = dpl::graph::vertex*;
  using edge_descriptor = dpl::graph::directed_edge*;

  using out_edge_iterator = dpl::graph::out_edge_iterator;

  using traversal_category = incidence_graph_tag;
  using directed_category = directed_tag;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = dpl::graph::dc_graph::vertices_size_type;
  using edges_size_type = dpl::graph::dc_graph::vertices_size_type;
  using degree_size_type = dpl::graph::dc_graph::degree_size_type;
};




