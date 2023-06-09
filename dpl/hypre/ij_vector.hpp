/*
 * This file is part of Rapid Reservoir Modelling Software.
 *   | https://rapidreservoir.org/
 *   | https://bitbucket.org/rapidreservoirmodelling/rrm/
 *
 * Copyright (c) 2022
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * RRM is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * RRM is distributed in the hope that it will be useful,                   
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

#include <tuple>

namespace dpl::hypre
{
#ifndef HYPRE_SEQUENTIAL
  inline std::tuple<HYPRE_BigInt, HYPRE_BigInt> mpi_part(HYPRE_BigInt nrows) {
    int w_size, w_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    
    auto count = nrows/w_size;
    auto remainder = nrows%w_size;
    HYPRE_BigInt start, stop;

    if (w_rank < remainder) {
      // The first 'remainder' ranks get 'count + 1' tasks each
      start = w_rank*(count + 1);
      stop = start + count;
    }
    else {
      // The remaining 'size - remainder' ranks get 'count' task each
      start = w_rank*count + remainder;
      stop = start + (count - 1);
    }

    return {start, stop};
  }
#endif
  

  
  class ij_vector
  {
    HYPRE_IJVector v_ = nullptr;

    void allocate_only(HYPRE_BigInt nvalues) {
    #ifdef HYPRE_SEQUENTIAL
      static constexpr auto comm = 0;
      static constexpr auto jlower = 0;
      const auto jupper = nvalues - 1;
    #else
      static constexpr auto comm = MPI_COMM_WORLD;
      auto [jlower, jupper] = mpi_part(nvalues);
    #endif

      HYPRE_IJVectorCreate(comm, jlower, jupper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);
      HYPRE_IJVectorAssemble(v_);
    }

    void allocate_assign(HYPRE_BigInt nvalues, const HYPRE_Complex* values) {
    #ifdef HYPRE_SEQUENTIAL
      static constexpr auto comm = 0;
      static constexpr auto jlower = 0;
      const auto jupper = nvalues - 1;
    #else
      static constexpr auto comm = MPI_COMM_WORLD;
      auto [jlower, jupper] = mpi_part(nvalues);
    #endif
      
      HYPRE_IJVectorCreate(comm, jlower, jupper, &v_);
      HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(v_);

      //      HYPRE_IJVectorSetValues(_v, nrows, rows, values);
      for (auto i = jlower; i <= jupper; ++i)
        HYPRE_IJVectorSetValues(v_, 1, &i, values + i);

      HYPRE_IJVectorAssemble(v_);
    }
  
  public:
    ij_vector(HYPRE_BigInt nvalues) {
      allocate_only(nvalues);
    }
  
    ij_vector(HYPRE_BigInt nvalues, const HYPRE_Complex* values) {
      allocate_assign(nvalues, values);
    }

    ~ij_vector() {
      if (v_)
        HYPRE_IJVectorDestroy(v_);
    }

    HYPRE_ParVector par_ref() const {
      HYPRE_ParVector pointer;
      HYPRE_IJVectorGetObject(v_, reinterpret_cast<void**>(&pointer));
      return pointer;
    }

    void get_values(HYPRE_Int nvalues, const HYPRE_BigInt* indices, HYPRE_Complex* values) const {
      HYPRE_IJVectorGetValues(v_, nvalues, indices, values);
    }

    ij_vector(const ij_vector& other) = delete;
    ij_vector(ij_vector&& other) = delete;
    ij_vector& operator=(const ij_vector& other) = delete;
    ij_vector& operator=(ij_vector&& other) = delete;
  };
}
