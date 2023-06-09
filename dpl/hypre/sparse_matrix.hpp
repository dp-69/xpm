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

#include <HYPRE_parcsr_ls.h>

#include <forward_list>
#include <vector>

namespace dpl::hypre
{
  class sparse_matrix
  {                    
    HYPRE_BigInt off_diag_count_;    
    
    friend struct InputDeprec;

  public:
    HYPRE_BigInt nrows;
    
    std::vector<HYPRE_Complex> diag;
    std::vector<std::forward_list<std::tuple<HYPRE_BigInt, HYPRE_Complex>>> off_diag;
    
    sparse_matrix() = default;
    
    explicit sparse_matrix(HYPRE_BigInt n) {
      nrows = n;
      diag.assign(n, 0);            
      // constants_.assign(n, 0);

      off_diag.resize(n);
      off_diag_count_ = 0;
    }

    void set_off_diag(HYPRE_BigInt i, HYPRE_BigInt j, HYPRE_Complex coef) {                    
      off_diag[i].emplace_front(j, coef);      
      ++off_diag_count_;
    }                  

    void add_diag(HYPRE_BigInt i, HYPRE_Complex value) {
      diag[i] += value;
    }

    auto& get_off_diag(HYPRE_BigInt i) {
      return off_diag[i];
    }

    void add_paired_diff_coef(HYPRE_BigInt i0, HYPRE_BigInt i1, HYPRE_Complex coef) {
      add_diag(i0, -coef);
      set_off_diag(i0, i1, coef);

      add_diag(i1, -coef);
      set_off_diag(i1, i0, coef);
    }       
  };
}
