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

#include <dpl/hypre/ij_matrix.hpp>

namespace dpl::hypre
{
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

    operator ls_known_ref() const {
      return { ncols.get(), cols.get(), values.get(), b.get() };
    }
  };

  class ls_known_storage_builder
  {
    HYPRE_BigInt nrows_;
    size_t nvalues_;
    ls_known_storage lks_;

    std::unique_ptr<size_t[]> diag_shift_;
    std::unique_ptr<HYPRE_Int[]> off_relative_;
  public:
    auto nvalues() const { return nvalues_; }

    void allocate_rows(HYPRE_BigInt nrows) {
      nrows_ = nrows;
      lks_.ncols = std::make_unique<HYPRE_Int[]>(nrows);
      lks_.b = std::make_unique<HYPRE_Complex[]>(nrows);

      for (HYPRE_BigInt i = 0; i < nrows; ++i) {
        lks_.ncols[i] = 1; // diag coef
        lks_.b[i] = 0;
      }

      nvalues_ = nrows;
    }

    void reserve_connection(HYPRE_BigInt i0, HYPRE_BigInt i1) {
      ++lks_.ncols[i0];
      ++lks_.ncols[i1];
      nvalues_ += 2;
    }

    void allocate_values() {
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

    void add_b(HYPRE_BigInt i, HYPRE_Complex value) {
      lks_.b[i] += value;
    }

    void add_diag(HYPRE_BigInt i, HYPRE_Complex value) {
      lks_.values[diag_shift_[i]] += value;
    }

    void set_connection(HYPRE_BigInt i0, HYPRE_BigInt i1, HYPRE_Complex value) {
      auto shift0 = diag_shift_[i0];
      lks_.values[shift0] += value; // diag
      shift0 += ++off_relative_[i0];
      lks_.cols[shift0] = i1;
      lks_.values[shift0] = -value;

      auto shift1 = diag_shift_[i1];
      lks_.values[shift1] += value; // diag
      shift1 += ++off_relative_[i1];
      lks_.cols[shift1] = i0;
      lks_.values[shift1] = -value;
    }

    auto acquire_storage() {
      return std::move(lks_);
    }
  };


  



  inline std::pair<HYPRE_Real, HYPRE_Int> solve(
    const index_range& range, const ls_known_ref& in, HYPRE_Complex* values,
    HYPRE_Real tolerance = 1.e-20, HYPRE_Int max_iterations = 20
  ) {

    HYPRE_Solver solver;
    
    HYPRE_BoomerAMGCreate(&solver);      

    HYPRE_BoomerAMGSetTol(solver, tolerance);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iterations);

    // HYPRE_BoomerAMGSetMaxLevels(solver, 50);
    // HYPRE_BoomerAMGSetPMaxElmts(solver, 0);
    // HYPRE_BoomerAMGSetMaxCoarseSize(solver, 18);

    // HYPRE_BoomerAMGSetCoarsenType(solver, 0);
    // HYPRE_BoomerAMGSetRestriction(solver, 2);
    // HYPRE_BoomerAMGSetTruncFactor(solver, 4);
    // HYPRE_BoomerAMGSetInterpType(solver, 6);
    // HYPRE_BoomerAMGSetStrongThreshold(solver, 0);
    // HYPRE_BoomerAMGSetTruncFactor(solver, 0);

    auto nrows = range.width();
    auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    std::iota(indices.get(), indices.get() + nrows, range.lower);

    ij_matrix A{range, indices.get(), in.ncols, in.cols, in.values};
    ij_vector b{range, in.b};
    ij_vector x{range};



    // const auto t0 = std::chrono::system_clock::now();

    // int w_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    //
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_0\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    HYPRE_BoomerAMGSetup(solver, A, b, x);


    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_1\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);
    
    HYPRE_BoomerAMGSolve(solver, A, b, x);
    

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_2\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);


    // const auto t1 = std::chrono::system_clock::now();
    // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";


    HYPRE_Real residual = 0;
    HYPRE_Int iters = 0;
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
    HYPRE_BoomerAMGGetNumIterations(solver, &iters);

    HYPRE_BoomerAMGDestroy(solver);                  
  
    x.get_values(nrows, indices.get(), values);

    return {residual, iters};
  }
}