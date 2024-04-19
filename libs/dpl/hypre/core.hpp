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
#include "../system.hpp"

#include <boost/iostreams/device/mapped_file.hpp>

namespace dpl::hypre
{
  inline constexpr auto hypre_input_name = "dpl-hypre-input";
  inline constexpr auto hypre_output_name = "dpl-hypre-output";

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

    // ls_known_storage() = default;
    // ls_known_storage(const ls_known_storage& other) = delete;
    // ls_known_storage(ls_known_storage&& other) noexcept = default;
    // ls_known_storage& operator=(const ls_known_storage& other) = delete;
    // ls_known_storage& operator=(ls_known_storage&& other) noexcept = default;
    // ~ls_known_storage() = default;

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
    std::size_t nvalues_;
    ls_known_storage lks_;

    std::unique_ptr<std::size_t[]> diag_shift_;
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

      diag_shift_ = std::make_unique<std::size_t[]>(nrows_);
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

  class ls_known_file_ref
  {
    class ptr_parser
    {
      void* ptr_;

    public:
      explicit ptr_parser(const char* ptr) : ptr_((void*) ptr) {}  // NOLINT(CppCStyleCast)

      template <typename T>
      auto& operator()(T& val) {
        val = *static_cast<T*>(ptr_);
        ptr_ = static_cast<char*>(ptr_) + sizeof(T);
        return *this;
      }
      
      template <typename T, typename S = int>
      auto& operator()(T*& val, S size = 0) {
        val = static_cast<T*>(ptr_);
        ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
        return *this;
      }

      auto* ptr() const { return ptr_; }
    };

    boost::iostreams::mapped_file_source mfs_;

  public:
    ls_known_file_ref() {   // NOLINT(cppcoreguidelines-pro-type-member-init)
      mfs_.open((std::filesystem::path("cache")/hypre_input_name).string());

      std::size_t nvalues;  // NOLINT(cppcoreguidelines-init-variables)

      ptr_parser{mfs_.data()}
        (nrows)
        (nvalues)
        (lkr.ncols, nrows)
        (lkr.b, nrows)
        (lkr.cols, nvalues)
        (lkr.values, nvalues)
        (tol)
        (max_iter)
        (agg_num_levels)
        (ranges);
    }

    ls_known_ref lkr;   
    HYPRE_BigInt nrows; 
    index_range* ranges;

    HYPRE_Real tol;
    HYPRE_Int max_iter;
    HYPRE_Int agg_num_levels;
  };

  using solve_result = std::tuple<std::unique_ptr<HYPRE_Complex[]>, HYPRE_Real, HYPRE_Int>;

  inline void save_input(
    ls_known_storage&& input, HYPRE_BigInt nrows, std::size_t nvalues, const std::vector<index_range>& blocks,  // NOLINT(cppcoreguidelines-rvalue-reference-param-not-moved)
    HYPRE_Real tol, HYPRE_Int max_iter, HYPRE_Int agg_num_levels)
  {
    stream_writer{std::filesystem::path{"cache"}/hypre_input_name}
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

    input.clear();
  }

  inline solve_result load_values(HYPRE_BigInt nrows) {
    auto values = std::make_unique<HYPRE_Complex[]>(nrows);
    HYPRE_Real residual;                                                  // NOLINT(cppcoreguidelines-init-variables)
    HYPRE_Int iters;                                                      // NOLINT(cppcoreguidelines-init-variables)

    stream_reader{std::filesystem::path{"cache"}/hypre_output_name}
      (values, nrows)
      (residual)
      (iters);

    return {std::move(values), residual, iters};
  }
  
  inline solve_result solve(std::unique_ptr<ls_known_file_ref> block, mpi::rank_t rank, HYPRE_Int print_level = 0)
  {
    #if (HYPRE_RELEASE_NUMBER == 22300)
      HYPRE_Init();
    #else
      HYPRE_Initialize();
    #endif

    HYPRE_Solver solver;          // NOLINT(cppcoreguidelines-init-variables)
    
    HYPRE_BoomerAMGCreate(&solver);      
    
    HYPRE_BoomerAMGSetTol(solver, block->tol);
    HYPRE_BoomerAMGSetMaxIter(solver, block->max_iter);
    HYPRE_BoomerAMGSetPrintLevel(solver, print_level);
    HYPRE_BoomerAMGSetAggNumLevels(solver, block->agg_num_levels);

    auto range = block->ranges[*rank];
    
    ij_matrix A{range, block->lkr.ncols, block->lkr.cols, block->lkr.values};
    ij_vector b{range, block->lkr.b};
    ij_vector x{range};

    block.reset();

    system::print_memory("<POST> HYPRE_Matrix Fill", print_level, rank);

    HYPRE_BoomerAMGSetup(solver, A, b, x);

    system::print_memory("<POST> HYPRE_BoomerAMGSetup", print_level, rank);

    HYPRE_BoomerAMGSolve(solver, A, b, x);

    system::print_memory("<POST> HYPRE_BoomerAMGSolve", print_level, rank);

    HYPRE_Real residual;          // NOLINT(cppcoreguidelines-init-variables)
    HYPRE_Int iters;              // NOLINT(cppcoreguidelines-init-variables)
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
    HYPRE_BoomerAMGGetNumIterations(solver, &iters);

    HYPRE_BoomerAMGDestroy(solver);                  

    auto values = std::make_unique<HYPRE_Complex[]>(range.width());
    x.get_values(range, values.get());

    HYPRE_Finalize();

    system::print_memory("<POST> HYPRE Destroy and Finalize", print_level, rank);

    return std::make_tuple(std::move(values), residual, iters);
  }


  inline solve_result solve(
    const ls_known_ref& in, const index_range& range, 
    HYPRE_Real tolerance = 1.e-20, HYPRE_Int max_iterations = 20/*, HYPRE_Int agg_num_levels = 0, HYPRE_Int print_level = 0*/)
  {
    #if (HYPRE_RELEASE_NUMBER == 22300)
      HYPRE_Init();
    #else
      HYPRE_Initialize();
    #endif

    HYPRE_Solver solver;          // NOLINT(cppcoreguidelines-init-variables)  
    
    HYPRE_BoomerAMGCreate(&solver);      

    HYPRE_BoomerAMGSetTol(solver, tolerance);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iterations);

    auto nrows = range.width();
    auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    std::iota(indices.get(), indices.get() + nrows, range.lower); // TODO: avoid indices allocation?

    ij_matrix A{range, indices.get(), in.ncols, in.cols, in.values};
    ij_vector b{range, in.b};
    ij_vector x{range};


    HYPRE_BoomerAMGSetup(solver, A, b, x);
    HYPRE_BoomerAMGSolve(solver, A, b, x);

    HYPRE_Real residual;          // NOLINT(cppcoreguidelines-init-variables)
    HYPRE_Int iters;              // NOLINT(cppcoreguidelines-init-variables)
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
    HYPRE_BoomerAMGGetNumIterations(solver, &iters);

    HYPRE_BoomerAMGDestroy(solver);                  

    auto values = std::make_unique<HYPRE_Complex[]>(range.width());
    // x.get_values(range, values.get());
    x.get_values(nrows, indices.get(), values.get());

    HYPRE_Finalize();

    return std::make_tuple(std::move(values), residual, iters);
  }

  inline void process() {
    mpi::rank_t rank;

    system::print_memory("START xpm MPI", 3, rank);

    auto input = std::make_unique<ls_known_file_ref>();
    auto global_offset = input->nrows*sizeof(HYPRE_Complex);
    auto range = input->ranges[*rank];

    auto [values, residual, iters] = solve(std::move(input), rank, 3);

    {
      MPI_File m_file;                                                              // NOLINT(cppcoreguidelines-init-variables)

      auto filename = (std::filesystem::path{"cache"}/hypre_output_name).string();

      MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &m_file);

      MPI_Status status;

      MPI_File_write_at_all(m_file, range.lower*sizeof(HYPRE_Complex),              // NOLINT(cppcoreguidelines-narrowing-conversions)
        values.get(), range.width(), 
        sizeof(HYPRE_Complex) == 8 ? MPI_DOUBLE : MPI_FLOAT,                        // NOLINT(CppUnreachableCode)
        &status);  

      if (rank) { /* root */
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
