/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2023
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

#include "ij_vector.hpp"

#include <numeric>

namespace dpl::hypre
{
  class ij_matrix
  {
    HYPRE_IJMatrix m_;

  public:
    ij_matrix(const index_range& range, HYPRE_Int* ncols, const HYPRE_BigInt* cols, const HYPRE_Complex* values) {  // NOLINT(cppcoreguidelines-pro-type-member-init)
      const auto [ilower, iupper] = range;
      const auto nrows = range.width();

      HYPRE_IJMatrixCreate(mpi::comm, ilower, iupper, ilower, iupper, &m_);
      HYPRE_IJMatrixSetObjectType(m_, HYPRE_PARCSR);

      auto shift = std::accumulate(ncols, ncols + ilower, static_cast<size_t>(0));

      {
        auto diag_sizes     = std::make_unique<HYPRE_Int[]>(nrows);
        auto off_diag_sizes = std::make_unique<HYPRE_Int[]>(nrows);
      
        std::fill_n(diag_sizes.get(), nrows, 0);
        std::fill_n(off_diag_sizes.get(), nrows, 0);

        const auto* cols_ptr = cols + shift;
        for (auto i = ilower; i <= iupper; ++i)
          for (HYPRE_Int col_idx = 0; col_idx < ncols[i]; ++col_idx)
            if (auto col = *cols_ptr++; ilower <= col && col <= iupper)
              ++diag_sizes[i - ilower];
            else
              ++off_diag_sizes[i - ilower];

        HYPRE_IJMatrixSetRowSizes(m_, ncols + ilower);
        HYPRE_IJMatrixSetDiagOffdSizes(m_, diag_sizes.get(), off_diag_sizes.get());
        HYPRE_IJMatrixSetMaxOffProcElmts(m_, 0);
      }

      HYPRE_IJMatrixInitialize(m_);
      for (auto i = ilower; i <= iupper; ++i) {
        HYPRE_IJMatrixSetValues(m_, 1, ncols + i, &i, cols + shift, values + shift);
        shift += ncols[i];
      }
      HYPRE_IJMatrixAssemble(m_);
    }

    ij_matrix(const index_range& range, const HYPRE_BigInt* rows, HYPRE_Int* ncols, const HYPRE_BigInt* cols, const HYPRE_Complex* values) {  // NOLINT(cppcoreguidelines-pro-type-member-init)
      const auto [ilower, iupper] = range;
      const auto nrows = range.width();

      HYPRE_IJMatrixCreate(mpi::comm, ilower, iupper, ilower, iupper, &m_);
      HYPRE_IJMatrixSetObjectType(m_, HYPRE_PARCSR);

      auto shift = std::accumulate(ncols, ncols + ilower, static_cast<size_t>(0));

      {
        auto diag_sizes =     std::make_unique<HYPRE_Int[]>(nrows);
        auto off_diag_sizes = std::make_unique<HYPRE_Int[]>(nrows);
      
        std::fill_n(diag_sizes.get(), nrows, 0);
        std::fill_n(off_diag_sizes.get(), nrows, 0);

        const auto* cols_ptr = cols + shift;
        for (auto i = ilower; i <= iupper; ++i)
          for (HYPRE_Int col_idx = 0; col_idx < ncols[i]; ++col_idx)
            if (auto col = *cols_ptr++; ilower <= col && col <= iupper)
              ++diag_sizes[i - ilower];
            else
              ++off_diag_sizes[i - ilower];

        HYPRE_IJMatrixSetRowSizes(m_, ncols + ilower);
        HYPRE_IJMatrixSetDiagOffdSizes(m_, diag_sizes.get(), off_diag_sizes.get());
        HYPRE_IJMatrixSetMaxOffProcElmts(m_, 0);
      }

      HYPRE_IJMatrixInitialize(m_);
      HYPRE_IJMatrixSetValues(m_, nrows, ncols + ilower, rows, cols + shift, values + shift);
      HYPRE_IJMatrixAssemble(m_);
    }

    ~ij_matrix() {
      if (m_)
        HYPRE_IJMatrixDestroy(m_);
    }

    ij_matrix(const ij_matrix& other) = delete;
    ij_matrix& operator=(const ij_matrix& other) = delete;
    ij_matrix(ij_matrix&& other) = delete;
    ij_matrix& operator=(ij_matrix&& other) = delete;

    operator HYPRE_ParCSRMatrix() const {             // NOLINT(CppNonExplicitConversionOperator)
      HYPRE_ParCSRMatrix ref;                         // NOLINT(cppcoreguidelines-init-variables)
      HYPRE_IJMatrixGetObject(m_, reinterpret_cast<void**>(&ref));
      return ref;
    }
  };
}