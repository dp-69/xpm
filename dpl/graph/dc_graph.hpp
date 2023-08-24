#pragma once

namespace dpl::graph::internal
{
  template<typename Pointer>
  class iterator_deference_pointer
  {
    Pointer ptr_;

  public:
    iterator_deference_pointer() = default;
    iterator_deference_pointer(Pointer ptr) : ptr_(ptr) {}

    using iterator_category = std::forward_iterator_tag;
    using value_type = Pointer;
    using difference_type = std::ptrdiff_t;
    // using distance_type = std::ptrdiff_t;
    // using pointer = Pointer;
    // using reference = Pointer;

    Pointer operator*() const { return ptr_; }
    // Pointer operator->() { return ptr; }

    iterator_deference_pointer& operator++() {
      ++ptr_;
      return *this;
    }  

    iterator_deference_pointer operator++(int) {
      iterator_deference_pointer temp = *this;
      ++*this;
      return temp;
    }

    friend bool operator==(const iterator_deference_pointer& l, const iterator_deference_pointer& r) { return l.ptr_ == r.ptr_; }
  };
}



namespace dpl::graph
{
  class vertex
  {
    friend class dc_properties;

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
    friend class dc_properties;

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


  using vertex_iterator = internal::iterator_deference_pointer<vertex*>;
  using out_edge_iterator = internal::iterator_deference_pointer<directed_edge*>;

  class dc_graph
  {
    std::size_t vertex_count_;

    std::unique_ptr<vertex[]> vertices_;
    std::unique_ptr<directed_edge[]> directed_edges_;
    std::unique_ptr<std::size_t[]> directed_edges_end_;

    // friend std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex*, const dc_graph&);
  public:
    using vertex_descriptor = vertex*;
    using edge_descriptor = directed_edge*;

    using vertices_size_type = std::size_t;
    using edges_size_type = std::size_t;
    using degree_size_type = std::size_t;

    explicit dc_graph(std::integral auto count)
      : vertex_count_(count),
        vertices_{std::make_unique<vertex[]>(count)},
        directed_edges_end_{std::make_unique<std::size_t[]>(count)} {}

    dc_graph(const dc_graph& other) = delete;
    dc_graph(dc_graph&& other) noexcept = delete;
    dc_graph& operator=(const dc_graph& other) = delete;
    dc_graph& operator=(dc_graph&& other) noexcept = delete;

    auto vertex_count() const { return vertex_count_; }

    const auto& vertices() const { return vertices_; }
    auto& vertices() { return vertices_; }
    vertex_descriptor get_vertex(std::integral auto idx) const { return vertices_.get() + idx; }

    const auto& directed_edges() const { return directed_edges_; }
    auto& directed_edges() { return directed_edges_; }
    edge_descriptor get_directed_edge(std::integral auto idx) const { return directed_edges_.get() + idx; }

    const auto& directed_edges_end() const { return directed_edges_end_; }
    auto& directed_edges_end() { return directed_edges_end_; }
  };

  inline dc_graph::vertices_size_type num_vertices(const dc_graph& g) {    
    return g.vertex_count();
  }

  inline vertex_iterator begin(dc_graph& g) { return g.vertices().get(); }
  inline vertex_iterator end(dc_graph& g) { return g.vertices().get() + g.vertex_count(); }

  inline std::pair<vertex_iterator, vertex_iterator> vertices(dc_graph& g) {
    return {begin(g), end(g)};
  }

  inline vertex* target(const directed_edge* e, const dc_graph&) {
    return e->v1;
  }

  inline void set_directed_edges_pair(vertex* v, vertex* u, directed_edge* vu, directed_edge* uv) {            
    vu->v1 = u;
    uv->v1 = v;

    vu->opposite = uv;
    uv->opposite = vu;
  }

  inline std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex* u, const dc_graph& g) {
    auto idx = u - g.vertices().get();

    return {
      g.get_directed_edge(idx == 0 ? 0 : g.directed_edges_end()[idx - 1]),
      g.get_directed_edge(g.directed_edges_end()[idx])
    };
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
  using vertex_descriptor = dpl::graph::dc_graph::vertex_descriptor;
  using edge_descriptor = dpl::graph::dc_graph::edge_descriptor;

  using out_edge_iterator = dpl::graph::out_edge_iterator;

  using traversal_category = incidence_graph_tag;
  using directed_category = directed_tag;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = dpl::graph::dc_graph::vertices_size_type;
  using edges_size_type = dpl::graph::dc_graph::vertices_size_type;
  using degree_size_type = dpl::graph::dc_graph::degree_size_type;
};




