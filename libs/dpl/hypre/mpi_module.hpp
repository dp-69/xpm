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

#include <dpl/hypre/core.hpp>

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace dpl::hypre::mpi
{
  class parser
  {
    void* ptr_;

  public:
    explicit parser(void* ptr) : ptr_(ptr) {}

    template <typename T>
    void read(T& val) {
      val = *static_cast<T*>(ptr_);
      ptr_ = static_cast<char*>(ptr_) + sizeof(T); 
    }
    
    template <typename T, typename S>
    void read_ref(T*& val, S size) {
      val = static_cast<T*>(ptr_);
      ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    }

    template <typename T, typename S>
    void read_copy(T* val, S size) {
      std::memcpy(val, ptr_, size*sizeof(T));
      ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    }

    template <typename T>
    void write(T val) {
      *static_cast<T*>(ptr_) = val;
      ptr_ = static_cast<char*>(ptr_) + sizeof(T);
    }

    template <typename T, typename S>
    void write(T* val, S size) {
      std::memcpy(ptr_, val, size*sizeof(T));
      ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    }

    auto* ptr() const { return ptr_; }
  };

  static inline std::filesystem::path mpi_exec;

  inline constexpr auto smo_hypre_input = "dpl-hypre-input";
  inline constexpr auto smo_hypre_output = "dpl-hypre-output";

  namespace bi = boost::interprocess;

  class block_info
  {
    bi::mapped_region region_;
    ls_known_ref lkr_;

  public:
    explicit block_info(int rank)
      : region_{bi::shared_memory_object{bi::open_only, smo_hypre_input, bi::read_only}, bi::read_only} {

      size_t nvalues;
      
      parser p{region_.get_address()};
      
      p.read(global_nrows);
      p.read(nvalues);
      p.read_ref(lkr_.ncols, global_nrows);
      p.read_ref(lkr_.b, global_nrows);
      p.read_ref(lkr_.cols, nvalues);
      p.read_ref(lkr_.values, nvalues);
      p.read(tol);
      p.read(max_iter);

      range = static_cast<index_range*>(p.ptr()) + rank;
    }

    HYPRE_BigInt global_nrows;
    index_range* range;

    HYPRE_Real tol;
    HYPRE_Int max_iter;

    operator ls_known_ref() const { return lkr_; }
  };

  using solve_result = std::tuple<std::unique_ptr<HYPRE_Complex[]>, HYPRE_Real, HYPRE_Int>;

  // struct solve_result
  // {
  //   std::unique_ptr<HYPRE_Complex[]> uptr;
  //   HYPRE_Real residual;
  //   HYPRE_Int iters;
  // };

  inline solve_result load_values(HYPRE_BigInt nrows) {
    auto values = std::make_unique<HYPRE_Complex[]>(nrows);
    HYPRE_Real residual;
    HYPRE_Int iters;

    auto region = bi::mapped_region{
      bi::shared_memory_object{bi::open_only, smo_hypre_output, bi::read_only},
      bi::read_only
    };

    parser p{region.get_address()};

    p.read_copy(values.get(), nrows);
    p.read(residual);
    p.read(iters);

    return {std::move(values), residual, iters};
  }

  inline void save_values(const HYPRE_Complex* values, HYPRE_BigInt nrows, HYPRE_Real residual, HYPRE_Int iters) {
    bi::shared_memory_object smo{bi::open_or_create, smo_hypre_output, bi::read_write};
    smo.truncate(nrows*sizeof(HYPRE_Complex) + sizeof(HYPRE_Real) + sizeof(HYPRE_Int));  // NOLINT(cppcoreguidelines-narrowing-conversions)

    auto region = bi::mapped_region{smo, bi::read_write};
    parser p{region.get_address()};
    p.write(values, nrows);
    p.write(residual);
    p.write(iters);
  }

  inline void save(
    const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<index_range>& blocks, HYPRE_Real tol, HYPRE_Int max_iter) {

    bi::shared_memory_object smo{bi::open_or_create, smo_hypre_input, bi::read_write};

    auto buffer_size = 
      sizeof(HYPRE_BigInt) + sizeof(size_t) +
      nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
      nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
      blocks.size()*sizeof(index_range) +
      sizeof(HYPRE_Real) + sizeof(HYPRE_Int)
    ;
    
    smo.truncate(buffer_size);  // NOLINT(cppcoreguidelines-narrowing-conversions)
    
    bi::mapped_region region(smo, bi::read_write);

    parser p{region.get_address()};
    p.write(nrows);
    p.write(nvalues);
    p.write(input.ncols, nrows);
    p.write(input.b, nrows);
    p.write(input.cols, nvalues);
    p.write(input.values, nvalues);
    p.write(tol);
    p.write(max_iter);
    
    std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(index_range));
  }

  inline void process() {
    int m_size, m_rank;  
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    static constexpr auto root = 0;

    auto input = block_info{m_rank};

    auto local_nrows = input.range->width();
    auto local_values = std::make_unique<HYPRE_Complex[]>(local_nrows);

    auto [residual, iters] = dpl::hypre::solve(*input.range, input, local_values.get(), input.tol, input.max_iter);

    std::unique_ptr<HYPRE_Complex[]> recvbuf;
    std::unique_ptr<int[]> recvcounts;
    std::unique_ptr<int[]> displs;

    if (m_rank == root) {
      recvbuf = std::make_unique<HYPRE_Complex[]>(input.global_nrows);
      recvcounts = std::make_unique<int[]>(m_size);
      displs = std::make_unique<int[]>(m_size);

      for (int i = 0; i < m_size; ++i) {
        recvcounts[i] = input.range[i].width();
        displs[i] = input.range[i].lower;
      }
    }

    static_assert(std::is_same_v<HYPRE_Complex, double>);

    MPI_Gatherv(
      local_values.get(), local_nrows, MPI_DOUBLE,
      recvbuf.get(), recvcounts.get(), displs.get(), MPI_DOUBLE,
      root, MPI_COMM_WORLD);

    if (m_rank == root)
      dpl::hypre::mpi::save_values(recvbuf.get(), input.global_nrows, residual, iters);
  }
}
