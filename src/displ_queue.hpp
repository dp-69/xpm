/*
 * This file is part of Extensive Pore Modelling (xpm).
 *   | https://github.com/dp-69/xpm
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * xpm is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * xpm is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with xpm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
    double r_cap;

    double pressure_cap() const {
      return 1/r_cap;
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
    void insert(macro_t i, double r_cap)  { set_.emplace(displ_elem::macro, *i, r_cap); }
    void insert(voxel_t i, double r_cap)  { set_.emplace(displ_elem::voxel, *i, r_cap); }
    void insert(throat_t i, double r_cap) { set_.emplace(displ_elem::throat, i, r_cap); }

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