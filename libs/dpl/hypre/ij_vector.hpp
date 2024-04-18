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

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

namespace dpl::hypre
{
  /**
   * \brief inclusive [lower, upper]
   */
  struct index_range
  {
    HYPRE_BigInt lower;
    HYPRE_BigInt upper;

    /**
     * \return
     *   local nrows (number of rows) as HYPRE_Int for hypre compatibility
     */
    HYPRE_Int width() const {
      return upper - lower + 1;  // NOLINT(cppcoreguidelines-narrowing-conversions)
    }
  };

  namespace mpi {
    #ifdef MPI_COMM_WORLD 
      inline MPI_Comm comm = MPI_COMM_WORLD;
    #else
      inline MPI_Comm comm = 0;
    #endif
  }
  
  class ij_vector
  {
    HYPRE_IJVector v_;
  
  public:
    explicit ij_vector(const index_range& range) {
      HYPRE_IJVectorCreate(mpi::comm, range.lower, range.upper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);
      HYPRE_IJVectorAssemble(v_);
    }

    explicit ij_vector(const index_range& range, const HYPRE_Complex* values) {
      const auto [jlower, jupper] = range;

      HYPRE_IJVectorCreate(mpi::comm, jlower, jupper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);

      for (auto i = jlower; i <= jupper; ++i)
        HYPRE_IJVectorSetValues(v_, 1, &i, values + i);

      HYPRE_IJVectorAssemble(v_);
    }

    ~ij_vector() {
      if (v_)
        HYPRE_IJVectorDestroy(v_);
    }

    ij_vector(const ij_vector& other) = delete;
    ij_vector(ij_vector&& other) = delete;
    ij_vector& operator=(const ij_vector& other) = delete;
    ij_vector& operator=(ij_vector&& other) = delete;

    operator HYPRE_ParVector() const {
      HYPRE_ParVector ref;
      HYPRE_IJVectorGetObject(v_, reinterpret_cast<void**>(&ref));
      return ref;
    }

    void get_values(const index_range& range, HYPRE_Complex* values) const {
      for (auto i = range.lower; i <= range.upper; ++i)
        HYPRE_IJVectorGetValues(v_, 1, &i, values++);
    }

    void get_values(HYPRE_Int nvalues, const HYPRE_BigInt* indices, HYPRE_Complex* values) const {
      HYPRE_IJVectorGetValues(v_, nvalues, indices, values);
    }
  };
}
