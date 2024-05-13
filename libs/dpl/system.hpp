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

#include "units.hpp"
#include "mpi_defs.hpp"

#include <fmt/core.h>

#include <iostream>

#ifdef _WIN32
  #include <windows.h>
  #undef max
  #undef min

  namespace dpl::system
  {
    template <typename T = units::byte>
    T get_memory_consumption() {
      MEMORYSTATUSEX status;
      status.dwLength = sizeof(status);

      GlobalMemoryStatusEx(&status);

      return units::byte{status.ullTotalPhys - status.ullAvailPhys};
    }
  }
#endif

namespace dpl::system
{
  inline void print_memory(const char* text, int print_level = 3, mpi::rank_t rank = mpi::rank_t::root()) {
    if (print_level > 0) {
      MPI_Barrier(mpi::comm);
      if (rank) { /* root */
        fmt::print("\n ~~~~~~~~~~~~~~~~ {} [{} GB] ~~~~~~~~~~~~~~~~\n", text, get_memory_consumption<units::gigabyte>());
        std::cout << "" << std::flush;
      }
      MPI_Barrier(mpi::comm);
    }
  }
}