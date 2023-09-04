#pragma once

#include <set>

#include "declarations.hpp"

namespace xpm
{
  struct occupancy_arrays
  {
    dpl::strong_array<net_tag, dpl::vector2d> macro;

    /**
     * \brief total, i.e. normal throat index
     */
    dpl::strong_array<size_t, double> throat; // TODO: if size_t, then should be simple array

    double occupancy(net_idx idx, double r_cap) {
      auto [c0, c1] = macro[idx];
      return c0 + c1*r_cap;
    }
  };

  enum struct displ_elem : unsigned char {
    macro,
    voxel,
    throat
  };

  struct displ_event
  {
    displ_elem elem;

    /**
     * \brief local, i.e. macro_idx or voxel_idx
     */
    size_t idx;
    double radius_cap;

    double pressure_cap() const {
      return 1/radius_cap;
    }
  };

  struct pressure_compare
  {
    bool operator()(const displ_event& l, const displ_event& r) const {
      return l.pressure_cap() < r.pressure_cap();  
    }
  };

  class displ_queue
  {
    std::multiset<displ_event, pressure_compare> set_;

  public:
    void insert(displ_elem elem, size_t idx, double r_cap) {
      set_.emplace(elem, idx, r_cap);
    }

    bool empty() const {
      return set_.empty();
    }

    auto& front() const {
      return *set_.begin();
    }

    void pop() {
      set_.erase(set_.begin());
    }
  };
}