#pragma once

#include <set>

#include "declarations.hpp"

namespace xpm
{
  enum struct displ_elem : unsigned char {
    macro,
    voxel,
    throat
  };

  struct displ_event
  {
    displ_elem elem;
    std::size_t local_idx; // local, i.e. macro_idx or voxel_idx
    double radius_cap;

    double pressure_cap() const {
      return 1/radius_cap;
    }
  };

  template<bool ascending>
  class displ_queue
  {
    struct comparator
    {
      bool operator()(const displ_event& l, const displ_event& r) const {
        return std::conditional_t<ascending, std::less<double>, std::greater<double>>{}(l.pressure_cap(), r.pressure_cap());
      }
    };

    std::multiset<displ_event, comparator> set_;

  public:
    void insert(macro_t i, double r_cap) { set_.emplace(displ_elem::macro, *i, r_cap); }
    void insert(voxel_t i, double r_cap) { set_.emplace(displ_elem::voxel, *i, r_cap); }
    void insert(std::size_t i, double r_cap) { set_.emplace(displ_elem::throat, i, r_cap); }

    bool empty() const {
      return set_.empty();
    }

    auto size() const {
      return set_.size();
    }

    auto& front() const {
      return *set_.begin();
    }

    void pop() {
      set_.erase(set_.begin());
    }
  };
}