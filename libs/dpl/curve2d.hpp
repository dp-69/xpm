/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "static_vector.hpp"

#include <algorithm>

namespace dpl
{
  /*
   * A sequence of 2D points
   * ordered by the first component
   * from low to high
   *
   */
  struct curve2d
  {
    std::vector<vector2d> storage;

    explicit operator bool() const {
      return !storage.empty();
    }

    auto inverse() const {
      if (storage.size() < 2)
        return *this;

      auto vec = storage;

      if (vec.front().y > vec.back().y)
        std::ranges::reverse(vec);

      vec.resize(
        std::ranges::unique(vec, {}, [](const vector2d& p) { return p.y; })
        .begin() - vec.begin());

      for (auto& p : vec)
        std::swap(p.x, p.y);
      return curve2d{std::move(vec)};
    }

    double solve(double value) const {
      return dpl::solve(std::span{storage}, value, flat);
    }

    double solve(double value, loga_t, linear_t) const {
      return dpl::solve(std::span{storage}, value, flat, loga, linear);
    }

    double operator()(double value) const {
      return solve(value);
    }

    double operator()(double value, loga_t, linear_t) const {
      return solve(value, loga, linear);
    }

    auto begin() const {
      return storage.begin();
    }

    auto end() const {
      return storage.end();
    }

    auto begin() {
      return storage.begin();
    }
    
    auto end() {
      return storage.end();
    }

    auto& front() const {
      return storage.front();
    }

    auto& back() const {
      return storage.back();
    }
  };

  // inline void unique(curve2d& curve) {
  //   auto& vec = curve.storage;
  //   vec.resize(std::ranges::unique(vec, {}, [](const vector2d& p) { return p.x; }).begin() - vec.begin());
  // }
}
