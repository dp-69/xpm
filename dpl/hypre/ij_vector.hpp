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
  namespace mpi
  {
    inline constexpr MPI_Comm comm =
      #ifdef MPI_COMM_WORLD 
        MPI_COMM_WORLD
      #else
        0
      #endif
    ;

    #ifndef HYPRE_SEQUENTIAL
    static std::pair<HYPRE_BigInt, HYPRE_BigInt> range; // TODO - needs to be reviwed
    #endif
  }
  
  class ij_vector
  {
    HYPRE_IJVector v_ = nullptr;
  
  public:
    ij_vector(HYPRE_BigInt jlower, HYPRE_BigInt jupper) {
      HYPRE_IJVectorCreate(mpi::comm, jlower, jupper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);
      HYPRE_IJVectorAssemble(v_);
    }

    ij_vector(HYPRE_BigInt jlower, HYPRE_BigInt jupper, const HYPRE_Complex* values) {
      HYPRE_IJVectorCreate(mpi::comm, jlower, jupper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);

      // auto count = jupper - jlower + 1;
      // auto indices = std::make_unique<HYPRE_BigInt[]>(count);
      // std::iota(indices.get(), indices.get() + count, jlower);
      // HYPRE_IJVectorSetValues(v_, count, indices.get(), values + jlower);
      // HYPRE_IJVectorSetValues(v_, jupper - jlower + 1, indices, values + jlower);

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

    HYPRE_ParVector par_ref() const {
      HYPRE_ParVector ref;  // NOLINT(cppcoreguidelines-init-variables)
      HYPRE_IJVectorGetObject(v_, reinterpret_cast<void**>(&ref));
      return ref;
    }

    void get_values(HYPRE_Int nvalues, const HYPRE_BigInt* indices, HYPRE_Complex* values) const {
      HYPRE_IJVectorGetValues(v_, nvalues, indices, values);
    }
  };
}
