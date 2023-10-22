#pragma once

#include "../general.hpp"

namespace dpl::graph
{
  class dc_graph
  {
    template<typename ValueType>
    class iterator
    {
      ValueType value_;

    public:
      iterator() = default;
      iterator(ValueType value) : value_(value) {}

      using iterator_category = std::forward_iterator_tag;
      using difference_type = std::ptrdiff_t;
      using value_type = ValueType;

      ValueType operator*() const {
        return value_;
      }

      auto& operator++() {
        ++value_;
        return *this;
      }  

      auto operator++(int) {
        auto temp = *this;
        ++*this;
        return temp;
      }

      friend bool operator==(const iterator& l, const iterator& r) { return l.value_ == r.value_; }
    };


  public:
    struct vertex_tag {};
    using vertex_t = strong_integer<std::int32_t, vertex_tag>;

    struct edge_tag {};
    using edge_t = strong_integer<std::size_t, edge_tag>;

  private:
    std::size_t vertex_count_;
    /**
     * \brief sparse ranges 
     */
    std::unique_ptr<edge_t[]> directed_edges_end_;

    std::size_t edge_count_;
    /**
     * \brief properties of the directed edges that are stored continuously in sparse format. directed_edges_end_ define partitioning
     */
    std::unique_ptr<edge_t[]> opposite_;
    std::unique_ptr<vertex_t[]> target_;

    /*
     * Euler tour related
     */
    std::unique_ptr<void*[]> vertex_entry_;
    std::unique_ptr<std::size_t[]> directed_edge_entry_;

    friend edge_t opposite(edge_t ab, const dc_graph& g);
    friend vertex_t target(edge_t ab, const dc_graph& g);

    friend class graph_generator_base;

  public:
    using vertex_iterator = iterator<vertex_t>;
    using out_edge_iterator = iterator<edge_t>;

    using vertices_size_type = std::size_t;
    using edges_size_type = std::size_t;
    using degree_size_type = std::size_t;

    dc_graph() = default;

    explicit dc_graph(std::integral auto vertex_count)  // NOLINT(cppcoreguidelines-pro-type-member-init)
      : vertex_count_(vertex_count),
        directed_edges_end_{std::make_unique<edge_t[]>(vertex_count)},
        vertex_entry_{std::make_unique<void*[]>(vertex_count)} {}

    dc_graph(const dc_graph& other) = delete;
    dc_graph(dc_graph&& other) noexcept = default;
    dc_graph& operator=(const dc_graph& other) = delete;
    dc_graph& operator=(dc_graph&& other) noexcept = default;

    void allocate_edges(std::size_t edge_count) {
      edge_count_ = edge_count;
      opposite_ = std::make_unique<edge_t[]>(edge_count*2);
      target_ = std::make_unique<vertex_t[]>(edge_count*2);
      directed_edge_entry_ = std::make_unique<std::size_t[]>(edge_count*2);
    }

    void*& vertex_entry(vertex_t i) const {
      return vertex_entry_[*i];
    }

    std::size_t& directed_edge_entry(edge_t i) const {
      return directed_edge_entry_[*i];
    }

    auto vertex_count() const {
      return vertex_count_;
    }

    auto edge_count() const {
      return edge_count_;
    }

    auto& directed_edges_end() {
      return directed_edges_end_;
    }

    out_edge_iterator out_edges_begin(vertex_t v) const {
      return v == 0 ? edge_t{0} : directed_edges_end_[*v - 1];
    }

    out_edge_iterator out_edges_end(vertex_t v) const {
      return directed_edges_end_[*v];
    }

    auto edges(vertex_t v) const {
      return std::ranges::subrange{out_edges_begin(v), out_edges_end(v)};
    }
  };

  inline dc_graph::vertices_size_type num_vertices(const dc_graph& g) {    
    return g.vertex_count();
  }

  inline std::pair<dc_graph::vertex_iterator, dc_graph::vertex_iterator> vertices(const dc_graph& g) {
    return {dc_graph::vertex_t{0}, dc_graph::vertex_t(g.vertex_count())};
  }

  inline dc_graph::vertex_t target(dc_graph::edge_t ab, const dc_graph& g) {
    return g.target_[*ab];
  }

  inline std::pair<dc_graph::out_edge_iterator, dc_graph::out_edge_iterator> out_edges(dc_graph::vertex_t v, const dc_graph& g) {
    return {g.out_edges_begin(v), g.out_edges_end(v)};
  }

  /**
   * \brief Boost
   */
  inline dc_graph::vertex_t source(dc_graph::edge_t, const dc_graph&) {}
  inline dc_graph::degree_size_type out_degree(dc_graph::vertex_t u, const dc_graph&) {}

  inline dc_graph::edge_t opposite(dc_graph::edge_t ab, const dc_graph& g) {
    return g.opposite_[*ab];
  }

  class graph_generator_base
  {
  protected:
    dc_graph g_;
    std::unique_ptr<std::size_t[]> shift_;
    std::size_t edge_count_ = 0;

    void set_dual(dc_graph::vertex_t a, dc_graph::vertex_t b, dc_graph::edge_t ab, dc_graph::edge_t ba) {            
      g_.target_[*ab] = b;
      g_.target_[*ba] = a;
      g_.opposite_[*ab] = ba;
      g_.opposite_[*ba] = ab;
    }

  public:
    explicit graph_generator_base(std::size_t vertex_count)
      : g_{vertex_count} {

      shift_= std::make_unique<std::size_t[]>(vertex_count + 1);
      std::fill_n(shift_.get(), vertex_count + 1, 0);
    }

    void allocate() {
      for (auto i = 0; i < g_.vertex_count(); ++i) {
        shift_[i + 1] += shift_[i];
        g_.directed_edges_end()[i] = dc_graph::edge_t{shift_[i + 1]};
      }

      g_.allocate_edges(edge_count_);
    }

    auto acquire() {
      return std::move(g_);
    }
  };

  template<typename Proj = std::identity>
  class graph_generator : public graph_generator_base
  {
    Proj proj_;

  public:
    explicit graph_generator(std::size_t vertex_count, Proj proj = {})
      : graph_generator_base(vertex_count), proj_{proj} {}

    void reserve(const auto l, const auto r) {
      ++shift_[*proj_(l) + 1]; 
      ++shift_[*proj_(r) + 1];
      edge_count_ += 1;
    }

    std::pair<dc_graph::edge_t, dc_graph::edge_t> set(const auto l, const auto r) {
      dc_graph::vertex_t proj_l = proj_(l);
      dc_graph::vertex_t proj_r = proj_(r);

      dc_graph::edge_t lr(shift_[*proj_l]);
      dc_graph::edge_t rl(shift_[*proj_r]);

      set_dual(proj_l, proj_r, lr, rl);

      ++shift_[*proj_l];
      ++shift_[*proj_r];

      return {lr, rl};
    }
  };
}


template <>
struct boost::graph_traits<dpl::graph::dc_graph>
{
  using vertex_descriptor = dpl::graph::dc_graph::vertex_t;
  using edge_descriptor = dpl::graph::dc_graph::edge_t;

  using out_edge_iterator = dpl::graph::dc_graph::out_edge_iterator;

  using traversal_category = incidence_graph_tag;
  using directed_category = undirected_tag;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = dpl::graph::dc_graph::vertices_size_type;
  using edges_size_type = dpl::graph::dc_graph::vertices_size_type;
  using degree_size_type = dpl::graph::dc_graph::degree_size_type;
};