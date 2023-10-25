#pragma once

#include "et_etnte_defs.hpp"
#include "dc_graph.hpp"

namespace dpl::graph
{
  class dc_traits
  {
    using et_ptr = et_traits::node_ptr;
    using et_cptr = et_traits::const_node_ptr;

    using vertex_t = dc_graph::vertex_t;
    using edge_t = dc_graph::edge_t;

    const dc_graph* g_;

  public:
    dc_traits() = default;

    explicit dc_traits(const dc_graph& graph)
      : g_(&graph) {}

    static edge_t get_directed_edge(et_cptr n) {
      return edge_t{mask_bit_balance::get_value<edge_t::value_type>(n->tag)};
    }

    static void set_directed_edge(et_ptr n, edge_t de) {
      mask_bit_balance::set_value(n->tag, *de);
    }

    static vertex_t get_vertex(et_cptr n) {
      return vertex_t{mask_bit_balance::get_value<vertex_t::value_type>(n->tag)};
    }

    static void set_vertex(et_ptr n, vertex_t v) {
      mask_bit_balance::set_value_bit(n->tag, *v);
    }                     

    /**
     * \brief checks if n is a vertex-type et entry
     */
    static bool is_loop_edge(et_cptr n) {
      return mask_bit_balance::get_bit(n->tag);
    }

    /**
     * \brief
     *   ordering of non-tree edges
     *   sorted by a pointing-in vertex, the pointing-out works as well
     */
    et_ptr get_ordering_vertex_entry(edge_t de) const {
      return get_entry(target(opposite(de, *g_), *g_));
    }

    et_ptr get_entry(vertex_t v) const {
      return static_cast<et_ptr>(g_->vertex_entry(v));
    }

    void set_entry(vertex_t v, et_ptr et) const {
      g_->vertex_entry(v) = et;
    }

    bool is_tree_edge(edge_t x) const {
      return !mask_bit::get_bit(g_->directed_edge_entry(x));
    }

    et_ptr get_tree_edge_entry(edge_t x) const {
      return mask_bit::get_ptr<et_ptr>(g_->directed_edge_entry(x));
    }

    void set_tree_edge_entry(edge_t x, et_ptr y) const {
      mask_bit::set_ptr(g_->directed_edge_entry(x), y);
    }

    void set_non_tree_edge(edge_t x) const { // TODO
      return mask_bit::set_bit(g_->directed_edge_entry(x));
    }

    bool is_null_entry(edge_t x) const {
      return g_->directed_edge_entry(x) == 0;
    }

    void set_null_entry(edge_t x) const {
      g_->directed_edge_entry(x) = 0;
    }
  };
}