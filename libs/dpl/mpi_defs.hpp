/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * DPL is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * DPL is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <mpi.h>

#include <filesystem>

namespace dpl::mpi
{
  static inline std::filesystem::path exec;

  #ifdef MPI_COMM_WORLD 
    inline MPI_Comm comm = MPI_COMM_WORLD;
  #else
    inline MPI_Comm comm = 0;
  #endif

  struct rank_t
  {
    int value;

    rank_t() {  // NOLINT(cppcoreguidelines-pro-type-member-init)
      MPI_Comm_rank(comm, &value);
    }

    explicit constexpr rank_t(int value)
      : value{value} {}

    constexpr auto operator*() const noexcept {
      return value;
    }

    explicit constexpr operator bool() const noexcept {
      return value == 0;
    }

    static constexpr rank_t root() noexcept {
      return rank_t{0};
    }
  };
}