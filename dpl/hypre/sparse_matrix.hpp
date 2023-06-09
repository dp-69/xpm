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
  class SparseMatrix
  {                    
    HYPRE_Int off_diag_count_;    
    
    friend struct Input;
    
  public:
    HYPRE_Int nrows;
    
    std::vector<std::forward_list<std::tuple<HYPRE_Int, HYPRE_Real>>> off_diag_coefs;
    std::vector<HYPRE_Real> diagonal;
    
    SparseMatrix() = default;
    
    explicit SparseMatrix(HYPRE_Int n) {
      SetSize(n);
    }

    void SetSize(HYPRE_Int n) {       
      nrows = n;
      diagonal.assign(n, 0);            
      // constants_.assign(n, 0);

      off_diag_coefs.resize(n);
      off_diag_count_ = 0;
    }

    void SetOffDiagCoef(HYPRE_Int row, HYPRE_Int col, double coef) {                    
      off_diag_coefs[row].emplace_front(col, coef);      
      ++off_diag_count_;
    }                  

    void AddDiagCoef(HYPRE_Int i, double value) {
      diagonal[i] += value;
    }

    auto& GetOffDiagCoefs(HYPRE_Int row) {
      return off_diag_coefs[row];
    }

    void AddDifferenceCoefs(HYPRE_Int owner, HYPRE_Int adj, HYPRE_Real coef) {
      AddDiagCoef(owner, -coef);
      SetOffDiagCoef(owner, adj, coef);

      AddDiagCoef(adj, -coef);
      SetOffDiagCoef(adj, owner, coef);
    }       
  };
}
