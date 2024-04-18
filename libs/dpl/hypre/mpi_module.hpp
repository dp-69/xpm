/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2023
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

#include "core.hpp"

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <filesystem>

namespace dpl::hypre::mpi
{
  

  class stream_writer
  {
    std::ofstream stream_;

  public:
    explicit stream_writer(const auto& path, std::ios_base::openmode mode = std::ios::binary/*, std::ofstream::pos_type pos = 0*/)
      : stream_{path, mode} {
      // stream_.seekp(pos);
    }

    template <typename T>
    auto& operator()(T val) {
      stream_.write(reinterpret_cast<const char*>(&val), sizeof(T));
      return *this;
    }

    template <typename T, typename S>
    auto& operator()(T* ptr, S size) {
      stream_.write(reinterpret_cast<const char*>(ptr), size*sizeof(T));
      return *this;
    }

    template <typename T, typename S>
    auto& operator()(const std::unique_ptr<T[]>& ptr, S size) {
      return (*this)(ptr.get(), size);
    }
  };

  static inline std::filesystem::path mpi_exec;
  static inline constexpr int root = 0;

  

  

  using solve_result = std::tuple<std::unique_ptr<HYPRE_Complex[]>, HYPRE_Real, HYPRE_Int>;

  inline solve_result load_values_file(HYPRE_BigInt nrows) {
    auto values = std::make_unique<HYPRE_Complex[]>(nrows);
    HYPRE_Real residual;
    HYPRE_Int iters;

    std::ifstream is{std::filesystem::path{"cache"}/smo_hypre_output, std::ios::binary};
    
    is.read(reinterpret_cast<char*>(values.get()), nrows*sizeof(HYPRE_Complex));  // NOLINT(cppcoreguidelines-narrowing-conversions)
    is.read(reinterpret_cast<char*>(&residual), sizeof(HYPRE_Real));
    is.read(reinterpret_cast<char*>(&iters), sizeof(HYPRE_Int));

    return {std::move(values), residual, iters};
  }

  // inline solve_result load_values(HYPRE_BigInt nrows) {
  //   auto values = std::make_unique<HYPRE_Complex[]>(nrows);
  //   HYPRE_Real residual;
  //   HYPRE_Int iters;
  //
  //   auto region = bi::mapped_region{
  //     bi::shared_memory_object{bi::open_only, smo_hypre_output, bi::read_only},
  //     bi::read_only
  //   };
  //
  //   parser p{region.get_address()};
  //
  //   p.read_copy(values.get(), nrows);
  //   p.read(residual);
  //   p.read(iters);
  //
  //   return {std::move(values), residual, iters};
  // }

  // inline void save_values(const HYPRE_Complex* values, HYPRE_BigInt nrows, HYPRE_Real residual, HYPRE_Int iters) {
  //   bi::shared_memory_object smo{bi::open_or_create, smo_hypre_output, bi::read_write};
  //   smo.truncate(nrows*sizeof(HYPRE_Complex) + sizeof(HYPRE_Real) + sizeof(HYPRE_Int));  // NOLINT(cppcoreguidelines-narrowing-conversions)
  //
  //   auto region = bi::mapped_region{smo, bi::read_write};
  //   parser p{region.get_address()};
  //   p.write(values, nrows);
  //   p.write(residual);
  //   p.write(iters);
  // }


  inline void save_and_reserve_file(
    ls_known_storage&& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<index_range>& blocks,  // NOLINT(cppcoreguidelines-rvalue-reference-param-not-moved)
    HYPRE_Real tol, HYPRE_Int max_iter, HYPRE_Int agg_num_levels)
  {
    stream_writer{std::filesystem::path{"cache"}/smo_hypre_input}
      (nrows)
      (nvalues)
      (input.ncols, nrows)
      (input.b, nrows)
      (input.cols, nvalues)
      (input.values, nvalues)
      (tol)
      (max_iter)
      (agg_num_levels)
      (blocks.data(), blocks.size());

    // {
    //   bi::shared_memory_object{bi::open_or_create, smo_hypre_output, bi::read_write}
    //     .truncate(nrows*sizeof(HYPRE_Complex) + sizeof(HYPRE_Real) + sizeof(HYPRE_Int));  // NOLINT(cppcoreguidelines-narrowing-conversions)
    // }

    input.clear();
  }

  // inline void save_and_reserve_file(
  //   const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<index_range>& blocks,
  //   HYPRE_Real tol, HYPRE_Int max_iter, HYPRE_Int agg_num_levels)
  // {
  //   parser_stream_writer{std::filesystem::path{"cache"}/smo_hypre_input}
  //     (nrows)
  //     (nvalues)
  //     (input.ncols, nrows)
  //     (input.b, nrows)
  //     (input.cols, nvalues)
  //     (input.values, nvalues)
  //     (tol)
  //     (max_iter)
  //     (agg_num_levels)
  //     (blocks.data(), blocks.size());
  //
  //   {
  //     bi::shared_memory_object{bi::open_or_create, smo_hypre_output, bi::read_write}
  //       .truncate(nrows*sizeof(HYPRE_Complex) + sizeof(HYPRE_Real) + sizeof(HYPRE_Int));  // NOLINT(cppcoreguidelines-narrowing-conversions)
  //   }
  // }

  // inline void save_and_reserve(
  //   const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<index_range>& blocks,
  //   HYPRE_Real tol, HYPRE_Int max_iter, HYPRE_Int agg_num_levels) {
  //
  //   bi::shared_memory_object smo{bi::open_or_create, smo_hypre_input, bi::read_write};
  //
  //   auto buffer_size = 
  //     sizeof(HYPRE_BigInt) + sizeof(size_t) +
  //     nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
  //     nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
  //     blocks.size()*sizeof(index_range) +
  //     sizeof(HYPRE_Real) + sizeof(HYPRE_Int) + sizeof(HYPRE_Int)
  //   ;
  //   
  //   smo.truncate(buffer_size);  // NOLINT(cppcoreguidelines-narrowing-conversions)
  //   
  //   bi::mapped_region region(smo, bi::read_write);
  //
  //   parser p{region.get_address()};
  //   p.write(nrows);
  //   p.write(nvalues);
  //   p.write(input.ncols, nrows);
  //   p.write(input.b, nrows);
  //   p.write(input.cols, nvalues);
  //   p.write(input.values, nvalues);
  //   p.write(tol);
  //   p.write(max_iter);
  //   p.write(agg_num_levels);
  //   
  //   std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(index_range));
  //
  //   {
  //     bi::shared_memory_object{bi::open_or_create, smo_hypre_output, bi::read_write}
  //       .truncate(nrows*sizeof(HYPRE_Complex) + sizeof(HYPRE_Real) + sizeof(HYPRE_Int));  // NOLINT(cppcoreguidelines-narrowing-conversions)
  //   }
  // }

  

  inline void process() {
    int /*m_size, */m_rank;                                                                       // NOLINT(cppcoreguidelines-init-variables)
    // MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    try_report_memory("START xpm MPI", 3, m_rank);

    auto input_uptr = std::make_unique<block_info>(m_rank);
    auto global_offset = input_uptr->global_nrows*sizeof(HYPRE_Complex);
    auto range = *input_uptr->range;

    auto [values, residual, iters] = solve(std::move(input_uptr), 3);

    {
      MPI_File m_file;                                                              // NOLINT(cppcoreguidelines-init-variables)

      auto filename = (std::filesystem::path{"cache"}/smo_hypre_output).string();

      MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &m_file);

      MPI_Status status;

      MPI_File_write_at_all(m_file, range.lower*sizeof(HYPRE_Complex),              // NOLINT(cppcoreguidelines-narrowing-conversions)
        values.get(), range.width(), 
        sizeof(HYPRE_Complex) == 8 ? MPI_DOUBLE : MPI_FLOAT,                        // NOLINT(CppUnreachableCode)
        &status);  

      if (m_rank == root) {
        MPI_File_write_at(m_file, global_offset,                      &residual, sizeof(HYPRE_Real), MPI_BYTE, &status);  // NOLINT(cppcoreguidelines-narrowing-conversions)
        MPI_File_write_at(m_file, global_offset + sizeof(HYPRE_Real), &iters,    sizeof(HYPRE_Int),  MPI_BYTE, &status);  // NOLINT(cppcoreguidelines-narrowing-conversions)
      }
      
      MPI_File_close(&m_file);

      MPI_Barrier(MPI_COMM_WORLD);

      // if (m_rank == root)
      //   stream_writer{filename, std::ios::binary | std::ios::app/*, global_count*sizeof(HYPRE_Complex)*/}
      //     (residual)
      //     (iters);
    }
  }
}
