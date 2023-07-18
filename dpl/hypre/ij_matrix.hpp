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

#include <dpl/hypre/ij_vector.hpp>

namespace dpl::hypre
{
  class ij_matrix
  {
    HYPRE_IJMatrix m_ = nullptr;

  public:
    void allocate_assign(HYPRE_BigInt nrows, HYPRE_Int* ncols, const HYPRE_BigInt* cols, const HYPRE_Complex* coefs) {
      #ifdef HYPRE_SEQUENTIAL
        static constexpr auto comm = 0;
        static constexpr auto ilower = 0;
        const auto iupper = nrows - 1;
      #else
        static constexpr auto comm = MPI_COMM_WORLD;
        auto [ilower, iupper] = dpl::hypre::mpi_block::range; // mpi_part(nrows);
      #endif
      
      HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &m_);
      
      HYPRE_IJMatrixSetObjectType(m_, HYPRE_PARCSR);
      // HYPRE_IJMatrixSetRowSizes(A, ncols.data());    

      // std::vector<HYPRE_Int> zeroArray(nrows, 0);
      // HYPRE_IJMatrixSetDiagOffdSizes(m_, ncols_per_row, zeroArray.data());

      HYPRE_IJMatrixInitialize(m_);
      // HYPRE_IJMatrixSetValues(_m, nrows, ncols, rows, cols, values);

      size_t shift = 0;
      
      for (HYPRE_BigInt i = 0; i < ilower; ++i)
        shift += ncols[i];
      
      for (HYPRE_BigInt i = ilower; i <= iupper; ++i) {
        HYPRE_IJMatrixSetValues(m_, 1, ncols + i, &i, cols + shift, coefs + shift);
        shift += ncols[i];
      }

      HYPRE_IJMatrixAssemble(m_);
    }

    auto par_ref() const {
      HYPRE_ParCSRMatrix ref;
      HYPRE_IJMatrixGetObject(m_, reinterpret_cast<void**>(&ref));
      return ref;
    }

    ij_matrix(HYPRE_BigInt nrows, HYPRE_Int* ncols, const HYPRE_BigInt* cols, const HYPRE_Complex* coefs) {
      allocate_assign(nrows, ncols, cols, coefs);
    }

    ~ij_matrix() {
      if (m_)
        HYPRE_IJMatrixDestroy(m_);
    }

    ij_matrix(const ij_matrix& other) = delete;
    ij_matrix& operator=(const ij_matrix& other) = delete;
    ij_matrix(ij_matrix&& other) = delete;
    ij_matrix& operator=(ij_matrix&& other) = delete;
  };
}