#pragma once

#include "../general.hpp"

namespace dpl::graph::internal
{
  template<typename Integral>
  class iterator_deference_integral
  {
    Integral value_;

  public:
    iterator_deference_integral() = default;
    iterator_deference_integral(Integral value) : value_(value) {}

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Integral;
    
    // using distance_type = std::ptrdiff_t;
    // using pointer = Pointer;
    // using reference = Pointer;

    Integral operator*() const { return value_; }
    // Pointer operator->() { return ptr; }

    auto& operator++() {
      ++value_;
      return *this;
    }  

    auto operator++(int) {
      auto temp = *this;
      ++*this;
      return temp;
    }

    friend bool operator==(const iterator_deference_integral& l, const iterator_deference_integral& r) { return l.value_ == r.value_; }
  };

  template<typename Ptr>
  class iterator_deference_pointer
  {
    Ptr ptr_;

  public:
    iterator_deference_pointer() = default;
    iterator_deference_pointer(Ptr ptr) : ptr_(ptr) {}

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Ptr;
    
    // using distance_type = std::ptrdiff_t;
    // using pointer = Pointer;
    // using reference = Pointer;

    Ptr operator*() const { return ptr_; }
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
  // class vertex
  // {
  //   friend class dc_properties;
  //
  // public:
  //   vertex() = default;
  //   vertex(const vertex& other) = delete;
  //   vertex(vertex&& other) noexcept = delete;
  //   vertex& operator=(const vertex& other) = delete;
  //   vertex& operator=(vertex&& other) noexcept = delete;
  // };

  using TEST_VERTEX_DESC_TYPE = int32_t;

  class directed_edge
  {
    friend class dc_properties;

    std::size_t entry_ = 0;


  public:
    directed_edge() = default;
    directed_edge(const directed_edge& other) = delete;
    directed_edge(directed_edge&& other) noexcept = delete;
    directed_edge& operator=(const directed_edge& other) = delete;
    directed_edge& operator=(directed_edge&& other) noexcept = delete;

    directed_edge* opposite;
    TEST_VERTEX_DESC_TYPE v1;  // pointing-to vertex  
  };


  using vertex_iterator = internal::iterator_deference_integral<TEST_VERTEX_DESC_TYPE>;//TEST_VERTEX_DESC_TYPE;  // TODO // internal::iterator_deference_pointer<vertex*>;
  using out_edge_iterator = internal::iterator_deference_pointer<directed_edge*>;

  class dc_graph
  {
    std::size_t vertex_count_;
    std::size_t edge_count_;

    std::unique_ptr<void*[]> vertex_entry_;

    std::unique_ptr<directed_edge[]> directed_edges_;
    std::unique_ptr<std::size_t[]> directed_edges_end_;

  public:
    using vertex_descriptor = TEST_VERTEX_DESC_TYPE;
    using edge_descriptor = directed_edge*;

    using vertices_size_type = std::size_t;
    using edges_size_type = std::size_t;
    using degree_size_type = std::size_t;

    dc_graph() = default;

    explicit dc_graph(std::integral auto vertex_count)
      : vertex_count_(vertex_count),
        vertex_entry_{std::make_unique<void*[]>(vertex_count_)},
        directed_edges_end_{std::make_unique<std::size_t[]>(vertex_count)} {}

    dc_graph(const dc_graph& other) = delete;
    dc_graph(dc_graph&& other) noexcept = default;
    dc_graph& operator=(const dc_graph& other) = delete;
    dc_graph& operator=(dc_graph&& other) noexcept = default;

    void allocate_edges(std::size_t edge_count) {
      edge_count_ = edge_count;
      directed_edges_ = std::make_unique<directed_edge[]>(edge_count*2);
    }

    auto vertex_count() const {
      return vertex_count_;
    }


    void* vertex_entry(vertex_descriptor i) const {
      return vertex_entry_[i];
    }

    void set_vertex_entry(vertex_descriptor i, void* entry) const {
      vertex_entry_[i] = entry;
    }

    vertex_iterator begin() const {
      return 0;
    }

    vertex_iterator end() const {
      return vertex_count_;
    }

    auto vertices() const {
      return std::ranges::subrange{begin(), end()};
    }

    // ---------------------------

    auto idx(const directed_edge* v) const {
      return v - directed_edges_.get();
    }

    auto edge_count() const {
      return edge_count_;
    }

    // const auto& directed_edges() const {
    //   return directed_edges_;
    // }

    edge_descriptor get_directed_edge(std::integral auto de_idx) const {
      return directed_edges_.get() + de_idx;
    }

    // const auto& directed_edges_end() const {
    //   return directed_edges_end_;
    // }

    auto& directed_edges_end() {
      return directed_edges_end_;
    }

    out_edge_iterator out_edges_begin(vertex_descriptor v) const {
      return &directed_edges_[v == 0 ? 0 : directed_edges_end_[v - 1]];
    }

    out_edge_iterator out_edges_end(vertex_descriptor v) const {
      return &directed_edges_[directed_edges_end_[v]];
    }

    auto edges(vertex_descriptor v) const {
      return std::ranges::subrange{out_edges_begin(v), out_edges_end(v)};
    }
  };

  inline dc_graph::vertices_size_type num_vertices(const dc_graph& g) {    
    return g.vertex_count();
  }

  

  inline std::pair<vertex_iterator, vertex_iterator> vertices(const dc_graph& g) {
    return {g.begin(), g.end()};
  }

  inline dc_graph::vertex_descriptor target(const directed_edge* e, const dc_graph&) {
    return e->v1;
  }

  inline std::pair<out_edge_iterator, out_edge_iterator> out_edges(const dc_graph::vertex_descriptor v, const dc_graph& g) {
    return {g.out_edges_begin(v), g.out_edges_end(v)};
  }
  
  inline dc_graph::vertex_descriptor source(const directed_edge*, const dc_graph&) {}
  inline dc_graph::degree_size_type out_degree(const dc_graph::vertex_descriptor u, const dc_graph&) {}


  






  template<typename Index = std::size_t>
  class graph_generator
  {
    dc_graph g_;
    std::unique_ptr<std::size_t[]> shift_;
    std::size_t edge_count_ = 0;

    using traits = strong_traits<Index>;

  public:
    explicit graph_generator(std::size_t vertex_count)
      : g_{vertex_count} {

      shift_= std::make_unique<std::size_t[]>(vertex_count + 1);
      std::fill_n(shift_.get(), vertex_count + 1, 0);
    }

    void reserve(Index l, Index r) {
      ++shift_[traits::get(l) + 1]; 
      ++shift_[traits::get(r) + 1];
      edge_count_ += 1;
    }

    void allocate() {
      for (auto i = 0; i < g_.vertex_count(); ++i) {
        shift_[i + 1] += shift_[i];
        g_.directed_edges_end()[i] = shift_[i + 1];
      }

      g_.allocate_edges(edge_count_);
    }

    void set(Index l, Index r) {
      auto set_pair = [](dc_graph::vertex_descriptor v, dc_graph::vertex_descriptor u, directed_edge* vu, directed_edge* uv) {            
        vu->v1 = u;
        uv->v1 = v;

        vu->opposite = uv;
        uv->opposite = vu;
      };

      set_pair(
        traits::get(l),
        traits::get(r),
        g_.get_directed_edge(shift_[traits::get(l)]++),
        g_.get_directed_edge(shift_[traits::get(r)]++));
    }

    auto acquire() {
      return std::move(g_);
    }
  };
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




