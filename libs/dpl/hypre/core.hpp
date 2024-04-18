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

#include "ij_matrix.hpp"
#include "../system.hpp" // TODO

#include <HYPRE_utilities.h>
#include <iostream>

#include <boost/iostreams/device/mapped_file.hpp>

#include <fmt/core.h> // TODO

namespace dpl::hypre
{
  inline constexpr auto smo_hypre_input = "dpl-hypre-input";
  inline constexpr auto smo_hypre_output = "dpl-hypre-output";



  inline void try_report_memory(const char* text, int print_level = 3, int m_rank = 0) {
    if (print_level > 0) {
      MPI_Barrier(mpi::comm);
      if (m_rank == 0) {
        fmt::print("\n ~~~~~~~~~~~~~~~~ {} [{} GB] ~~~~~~~~~~~~~~~~\n", text, get_memory_consumption<units::gigabyte>());
        std::cout << "" << std::flush;
      }
      MPI_Barrier(mpi::comm);
    }
  }






  /**
   * \brief
   *    nrows = length(ncols) = length(b), number of rows
   *    nvalues = sum(ncols) = length(cols) = length(values), number of non-zero coefficients
   */
  struct ls_known_ref
  {
    HYPRE_Int* ncols;       // number of non-zero columns
    HYPRE_BigInt* cols;     // non-zero columns
    HYPRE_Complex* values;  // non-zero coefs

    HYPRE_Complex* b;       // constant terms
  };

  struct ls_known_storage
  {
    std::unique_ptr<HYPRE_Int[]> ncols;       
    std::unique_ptr<HYPRE_BigInt[]> cols;     
    std::unique_ptr<HYPRE_Complex[]> values;  

    std::unique_ptr<HYPRE_Complex[]> b;

    ls_known_storage() = default;
    ls_known_storage(const ls_known_storage& other) = delete;
    ls_known_storage(ls_known_storage&& other) noexcept = default;
    ls_known_storage& operator=(const ls_known_storage& other) = delete;
    ls_known_storage& operator=(ls_known_storage&& other) noexcept = default;
    ~ls_known_storage() = default;

    operator ls_known_ref() const { // NOLINT(CppNonExplicitConversionOperator)
      return { ncols.get(), cols.get(), values.get(), b.get() };
    }

    void clear() {
      ncols.reset();
      cols.reset();
      values.reset();
      b.reset();
    }
  };

  template <typename Proj = std::identity>
  class ls_known_storage_builder
  {
    Proj proj_;

    HYPRE_BigInt nrows_;
    size_t nvalues_;
    ls_known_storage lks_;

    std::unique_ptr<size_t[]> diag_shift_;
    std::unique_ptr<HYPRE_Int[]> off_relative_;

  public:
    explicit ls_known_storage_builder(HYPRE_BigInt nrows, Proj proj = {})
      : proj_(proj) {
      nrows_ = nrows;
      lks_.ncols = std::make_unique<HYPRE_Int[]>(nrows);
      lks_.b = std::make_unique<HYPRE_Complex[]>(nrows);

      for (HYPRE_BigInt i = 0; i < nrows; ++i) {
        lks_.ncols[i] = 1; // diag coef
        lks_.b[i] = 0;
      }

      nvalues_ = nrows;
    }

    void reserve(/*HYPRE_BigInt*/auto i0, /*HYPRE_BigInt*/auto i1) {
      ++lks_.ncols[proj_(i0)];
      ++lks_.ncols[proj_(i1)];
      nvalues_ += 2;
    }

    void allocate() {
      lks_.cols = std::make_unique<HYPRE_BigInt[]>(nvalues_);
      lks_.values = std::make_unique<HYPRE_Complex[]>(nvalues_);

      diag_shift_ = std::make_unique<size_t[]>(nrows_);
      off_relative_ = std::make_unique<HYPRE_Int[]>(nrows_);

      diag_shift_[0] = 0;
      off_relative_[0] = 0; // pre increment will be required
      lks_.cols[0] = 0;
      lks_.values[0] = 0; // diag
      
      for (HYPRE_BigInt i = 1; i < nrows_; ++i) {
        auto s = diag_shift_[i - 1] + lks_.ncols[i - 1];
        diag_shift_[i] = s;
        off_relative_[i] = 0; // pre increment will be required
        lks_.cols[s] = i;
        lks_.values[s] = 0; // diag
      }
    }

    void add_b(/*HYPRE_BigInt*/auto i, HYPRE_Complex value) {
      lks_.b[proj_(i)] += value;
    }

    void add_diag(/*HYPRE_BigInt*/auto i, HYPRE_Complex value) {
      lks_.values[diag_shift_[proj_(i)]] += value;
    }

    void set(/*HYPRE_BigInt*/auto i0, /*HYPRE_BigInt*/auto i1, HYPRE_Complex value) {
      auto shift0 = diag_shift_[proj_(i0)];
      lks_.values[shift0] += value; // diag
      shift0 += ++off_relative_[proj_(i0)];
      lks_.cols[shift0] = proj_(i1);
      lks_.values[shift0] = -value;

      auto shift1 = diag_shift_[proj_(i1)];
      lks_.values[shift1] += value; // diag
      shift1 += ++off_relative_[proj_(i1)];
      lks_.cols[shift1] = proj_(i0);
      lks_.values[shift1] = -value;
    }

    auto acquire() {
      return std::move(lks_);
    }

    auto nvalues() const {
      return nvalues_;
    }
  };


  class parser
  {
    void* ptr_;

  public:
    explicit parser(void* ptr, std::size_t offset = 0) : ptr_(static_cast<unsigned char*>(ptr) + offset) {}

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

    // template <typename T, typename S>
    // void read_copy(T* val, S size) {
    //   std::memcpy(val, ptr_, size*sizeof(T));
    //   ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    // }

    // template <typename T>
    // void write(T val) {
    //   *static_cast<T*>(ptr_) = val;
    //   ptr_ = static_cast<char*>(ptr_) + sizeof(T);
    // }

    // template <typename T, typename S>
    // void write(T* val, S size) {
    //   std::memcpy(ptr_, val, size*sizeof(T));
    //   ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    // }

    auto* ptr() const { return ptr_; }

    // void advance(auto offset) {
    //   ptr_ = static_cast<char*>(ptr_) + offset;
    // }
  };

  class block_info
  {
    boost::iostreams::mapped_file_source mfs_;

    

  public:


    explicit block_info(int rank) {
      mfs_.open((std::filesystem::path("cache")/smo_hypre_input).string());
      // mfs = boost::iostreams::mapped_file_source{filename};

      size_t nvalues;

      parser p{(void*)mfs_.data()};  // NOLINT(CppCStyleCast)
      
      p.read(global_nrows);
      p.read(nvalues);
      p.read_ref(lkr.ncols, global_nrows);
      p.read_ref(lkr.b, global_nrows);
      p.read_ref(lkr.cols, nvalues);
      p.read_ref(lkr.values, nvalues);
      p.read(tol);
      p.read(max_iter);
      p.read(agg_num_levels);

      range = static_cast<index_range*>(p.ptr()) + rank;
    }

    ls_known_ref lkr;

    HYPRE_BigInt global_nrows;
    index_range* range;

    HYPRE_Real tol;
    HYPRE_Int max_iter;
    HYPRE_Int agg_num_levels;

    operator ls_known_ref() const { return lkr; }  // NOLINT(CppNonExplicitConversionOperator)
  };


  inline std::tuple<
    std::unique_ptr<HYPRE_Complex[]>,
    HYPRE_Real, HYPRE_Int>
  solve(std::unique_ptr<block_info> block, HYPRE_Int print_level = 0)
  {
    #if (HYPRE_RELEASE_NUMBER == 22300)
      HYPRE_Init();
    #else
      HYPRE_Initialize();
    #endif

    HYPRE_Solver solver;
    
    HYPRE_BoomerAMGCreate(&solver);      
    
    HYPRE_BoomerAMGSetTol(solver, block->tol);
    HYPRE_BoomerAMGSetMaxIter(solver, block->max_iter);
    HYPRE_BoomerAMGSetPrintLevel(solver, print_level);
    HYPRE_BoomerAMGSetAggNumLevels(solver, block->agg_num_levels);

    auto range = *block->range;
    // auto nrows = range.width();
    // auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    // std::iota(indices.get(), indices.get() + nrows, range.lower); // TODO: avoid indices allocation?

    
    ij_matrix A{range, /*indices.get(), */block->lkr.ncols, block->lkr.cols, block->lkr.values};
    ij_vector b{range, block->lkr.b};
    ij_vector x{range};

    block.reset();

    int m_rank;  
    MPI_Comm_rank(mpi::comm, &m_rank);

    try_report_memory("<POST> HYPRE_Matrix Fill", print_level, m_rank);

    HYPRE_BoomerAMGSetup(solver, A, b, x);

    try_report_memory("<POST> HYPRE_BoomerAMGSetup", print_level, m_rank);

    HYPRE_BoomerAMGSolve(solver, A, b, x);

    try_report_memory("<POST> HYPRE_BoomerAMGSolve", print_level, m_rank);

    HYPRE_Real residual = 0;
    HYPRE_Int iters = 0;
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
    HYPRE_BoomerAMGGetNumIterations(solver, &iters);

    HYPRE_BoomerAMGDestroy(solver);                  

    auto values = std::make_unique<HYPRE_Complex[]>(range.width());
    x.get_values(range, values.get());
    // x.get_values(nrows, indices.get(), values.get());

    HYPRE_Finalize();

    try_report_memory("<POST> HYPRE Destroy and Finalize", print_level, m_rank);

    return std::make_tuple(std::move(values), residual, iters);
  }


  inline std::tuple<
    std::unique_ptr<HYPRE_Complex[]>,
    HYPRE_Real, HYPRE_Int>
  solve(
    const ls_known_ref& in, const index_range& range, 
    HYPRE_Real tolerance = 1.e-20, HYPRE_Int max_iterations = 20, HYPRE_Int agg_num_levels = 0, HYPRE_Int print_level = 0)
  {
    #if (HYPRE_RELEASE_NUMBER == 22300)
      HYPRE_Init();
    #else
      HYPRE_Initialize();
    #endif

    HYPRE_Solver solver;
    
    HYPRE_BoomerAMGCreate(&solver);      

    HYPRE_BoomerAMGSetTol(solver, tolerance);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iterations);
    HYPRE_BoomerAMGSetPrintLevel(solver, print_level);

    HYPRE_BoomerAMGSetAggNumLevels(solver, agg_num_levels);
    // HYPRE_BoomerAMGSetAggInterpType(solver, 1); // ?

    // HYPRE_BoomerAMGSetStrongThreshold(solver, 0.1);

    // HYPRE_BoomerAMGSetRelaxType(solver, 9); 
    // HYPRE_BoomerAMGSetRelaxType()
    // HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    // HYPRE_BoomerAMGSetMaxLevels(solver, 3);

    auto nrows = range.width();
    auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    std::iota(indices.get(), indices.get() + nrows, range.lower); // TODO: avoid indices allocation?

    ij_matrix A{range, indices.get(), in.ncols, in.cols, in.values};
    ij_vector b{range, in.b};
    ij_vector x{range};


    int m_rank;  
    MPI_Comm_rank(mpi::comm, &m_rank);

    try_report_memory("<POST> HYPRE_Matrix Fill", print_level, m_rank);

    HYPRE_BoomerAMGSetup(solver, A, b, x);

    try_report_memory("<POST> HYPRE_BoomerAMGSetup", print_level, m_rank);

    HYPRE_BoomerAMGSolve(solver, A, b, x);

    try_report_memory("<POST> HYPRE_BoomerAMGSolve", print_level, m_rank);

    HYPRE_Real residual = 0;
    HYPRE_Int iters = 0;
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
    HYPRE_BoomerAMGGetNumIterations(solver, &iters);

    HYPRE_BoomerAMGDestroy(solver);                  

    auto values = std::make_unique<HYPRE_Complex[]>(range.width());
    // x.get_values(range, values.get());
    x.get_values(nrows, indices.get(), values.get());

    HYPRE_Finalize();

    try_report_memory("<POST> HYPRE Destroy and Finalize", print_level, m_rank);

    return std::make_tuple(std::move(values), residual, iters);
  }
}
