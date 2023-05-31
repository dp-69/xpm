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

#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif

// #include "rrm/fd/core/linear_solver/HypreVectorMatrix.hpp"

#include <boost/format/format_fwd.hpp>
#include <boost/format.hpp>

#include <vector>
#include <fstream>

#include <forward_list>











namespace dpl::hypre
{
  namespace internal
  {
    class Vector
    {
      HYPRE_IJVector v_ = nullptr;

    public:
      // HypreVector() = default;

      Vector(HYPRE_Int nrows/*, const HYPRE_Int* rows*/, const HYPRE_Real* values) {
        SetValues(nrows/*, rows*/, values);
      }

      Vector(HYPRE_Int nrows) {
        AllocateSpace(nrows);
      }

      ~Vector() {
        if (v_)
          HYPRE_IJVectorDestroy(v_);
      }

      void AllocateSpace(HYPRE_Int nrows) {
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, nrows - 1, &v_);
        HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(v_);

        HYPRE_IJVectorAssemble(v_);
      }

      void SetValues(HYPRE_Int nrows/*, const HYPRE_Int* rows*/, const HYPRE_Real* values) {
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, nrows - 1, &v_);
        HYPRE_IJVectorSetObjectType(v_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(v_);
        
        //      HYPRE_IJVectorSetValues(_v, nrows, rows, values);
        for (HYPRE_Int i = 0; i < nrows; ++i)
          HYPRE_IJVectorSetValues(v_, 1, &i, values + i);

        HYPRE_IJVectorAssemble(v_);
      }

      HYPRE_ParVector GetParVector() const {
        HYPRE_ParVector pointer;
        HYPRE_IJVectorGetObject(v_, reinterpret_cast<void**>(&pointer));
        return pointer;
      }

      void GetValues(HYPRE_Int nrows, const HYPRE_Int* rows, HYPRE_Real* values) const {
        HYPRE_IJVectorGetValues(v_, nrows, rows, values);
      }


      Vector(const Vector& other) = delete;
      Vector(Vector&& other) = delete;
      Vector& operator=(const Vector& other) = delete;
      Vector& operator=(Vector&& other) = delete;
    };

    class Matrix
    {
      HYPRE_IJMatrix m_ = nullptr;

    public:
      void SetValues(HYPRE_Int nrows, HYPRE_Int ncols, HYPRE_Int* ncols_per_row, const HYPRE_Int* cols_of_coefs, const HYPRE_Real* coefs) {
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, nrows - 1, 0, ncols - 1, &m_);

        HYPRE_IJMatrixSetObjectType(m_, HYPRE_PARCSR);
        // HYPRE_IJMatrixSetRowSizes(A, ncols.data());    

        // std::vector<HYPRE_Int> zeroArray(nrows, 0);
        // HYPRE_IJMatrixSetDiagOffdSizes(m_, ncols_per_row, zeroArray.data());

        HYPRE_IJMatrixInitialize(m_);
        // HYPRE_IJMatrixSetValues(_m, nrows, ncols, rows, cols, values);

        for (HYPRE_Int i = 0, shift = 0; i < nrows; ++i) {
          HYPRE_IJMatrixSetValues(m_, 1, ncols_per_row + i, &i, cols_of_coefs + shift, coefs + shift);
          shift += ncols_per_row[i];
        }

        HYPRE_IJMatrixAssemble(m_);
      }

      HYPRE_ParCSRMatrix GetParCSRMatrix() const {
        HYPRE_ParCSRMatrix csr;
        HYPRE_IJMatrixGetObject(m_, reinterpret_cast<void**>(&csr));
        return csr;
      }

      Matrix(HYPRE_Int nrows, HYPRE_Int ncols, HYPRE_Int* ncols_per_row, const HYPRE_Int* cols_of_coefs, const HYPRE_Real* coefs) {
        SetValues(nrows, ncols, ncols_per_row, cols_of_coefs, coefs);
      }

      // HypreMatrix() = default;

      ~Matrix() {
        if (m_)
          HYPRE_IJMatrixDestroy(m_);
      }

      Matrix(const Matrix& other) = delete;
      Matrix& operator=(const Matrix& other) = delete;
      Matrix(Matrix&& other) = delete;
      Matrix& operator=(Matrix&& other) = delete;
    };
  }













  
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
  
  struct Input
  {
    HYPRE_Int nrows;    
    std::vector<HYPRE_Int> ncols_per_row; // size of nrows       
    std::vector<HYPRE_Real> constants; // size of nrows

    std::vector<HYPRE_Int> cols_of_coefs; 
    std::vector<HYPRE_Real> coefs;



    Input() = default;

    Input(const Input& other) = default;
    Input(Input&& other) noexcept = default;
    Input& operator=(const Input& other) = default;
    Input& operator=(Input&& other) noexcept = default;
    

    Input(SparseMatrix& m, std::vector<double>&& free_terms) {             
      nrows = m.nrows;
      ncols_per_row.assign(m.nrows, 1);
      constants = std::move(free_terms);            
      
      const auto coef_count = nrows + m.off_diag_count_;
      
      cols_of_coefs.resize(coef_count);
      coefs.assign(coef_count, 0);

            
      for (HYPRE_Int row_idx = 0, coef_idx = 0; row_idx < nrows; row_idx++) {                
        auto diag_idx = coef_idx++;
        cols_of_coefs[diag_idx] = row_idx;
        coefs[diag_idx] = m.diagonal[row_idx];        
        
        for (auto [j, j_off_coef] : m.off_diag_coefs[row_idx]) {
          ++ncols_per_row[row_idx];
          
          auto off_diag_idx = coef_idx++;
          cols_of_coefs[off_diag_idx] = j;          
          coefs[off_diag_idx] = j_off_coef;          
        }                                                
      }          
    }
    
    void Export(const std::string& filename, bool numerical = false) {
      std::ofstream out(filename);      

      auto* colPtr = cols_of_coefs.data();
      auto* coefPtr = coefs.data();
      for (auto i = 0; i < nrows; ++i) {
        std::vector<HYPRE_Real> row_values(nrows, 0);

        for (auto j = 0; j < ncols_per_row[i]; ++j)
          row_values[*colPtr++] = *coefPtr++;

        for (auto j = 0; j < nrows; ++j) {          
          if (numerical)
            out << boost::format("%10.2e") % row_values[j];
          else
            out << (row_values[j] ? row_values[j] == 1 ? '1' : 'v' : '0'); // NOLINT(clang-diagnostic-float-conversion, clang-diagnostic-float-equal)
          
          out << ' ';
        }

        out << "| ";
                
        if (numerical)
          out << boost::format("%10.2e") % constants[i];
        else
          out << (constants[i] ? constants[i] == 1 ? '1' : 'v' : '0');  // NOLINT(clang-diagnostic-float-conversion, clang-diagnostic-float-equal)
        
        out << std::endl;
      }
      
      // out << "one" << " two" << boost::format("%e %.3e") % 1.5e3 % 0.0048799;

      out.close();
    }
    
    void Solve(HYPRE_Int* indices_ptr, HYPRE_Int indices_count, HYPRE_Real* output_ptr) {                  
      HYPRE_Solver solver;

      
      HYPRE_BoomerAMGCreate(&solver);      

      // (Optional) Set the convergence tolerance, if BoomerAMG is used as a solver. If it is used as a preconditioner,
      //   it should be set to 0. The default is 1.e-7.

      // HYPRE_Real tolerance = 1e-20;
      HYPRE_BoomerAMGSetTol(solver, 1.e-20);

      
      // (Optional) Sets AMG strength threshold. The default is 0.25. For 2d Laplace operators, 0.25 is a good value,
      //   for 3d Laplace operators, 0.5 or 0.6 is a better value. For elasticity problems, a large strength threshold,
      //     such as 0.9, is often better.
      // HYPRE_BoomerAMGSetStrongThreshold(solver, 0.5); // 0.7

      // 0 CLJP-coarsening (a parallel coarsening algorithm using independent sets.
      //   1 classical Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
      //   3 classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse
      //   points on the boundaries
      //   6 Falgout coarsening (uses 1 first, followed by CLJP using the interior coarse points
      //     generated by 1 as its first independent set)
      //   7 CLJP-coarsening (using a fixed random vector, for debugging purposes only)
      //   8 PMIS-coarsening (a parallel coarsening algorithm using independent sets, generating
      //     lower complexities than CLJP, might also lead to slower convergence)
      //   9 PMIS-coarsening (using a fixed random vector, for debugging purposes only)
      //   10 HMIS-coarsening (uses one pass Ruge-Stueben on each processor independently, followed
      //     by PMIS using the interior C-points generated as its first independent set)
      //   11 one-pass Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
      //   21 CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer
      //   22 CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer

      // The default is 6.

      
      // HYPRE_BoomerAMGSetCoarsenType(solver, 9); //HMIS - 10 PMIS - 9



      // (Optional) Defines which parallel interpolation operator is used. There are the following options for interp type:

       // 0 classical modified interpolation
       //  1 LS interpolation (for use with GSMG)
       //  2 classical modified interpolation for hyperbolic PDEs
       //  3 direct interpolation (with separation of weights)
       //  4 multipass interpolation
       //  5 multipass interpolation (with separation of weights)
       //  6 extended+i interpolation
       //  7 extended+i (if no common C neighbor) interpolation
       //  8 standard interpolation
       //  9 standard interpolation (with separation of weights)
       //  10 classical block interpolation (for use with nodal systems version only)
       //  11 classical block interpolation (for use with nodal systems version only)
       //  with diagonalized diagonal blocks
       //  12 FF interpolation
       //  13 FF1 interpolation
       //  14 extended interpolation

      // The default is 0.

      // HYPRE_BoomerAMGSetInterpType(solver, 8);

  

      
      const internal::Matrix A_matrix(nrows, nrows, ncols_per_row.data(), cols_of_coefs.data(), coefs.data());
      const internal::Vector B_vector(nrows, constants.data());
      const internal::Vector x_vector(nrows);
      
      auto* A_parcsr = A_matrix.GetParCSRMatrix();
      auto* B_par = B_vector.GetParVector();
      auto* x_par = x_vector.GetParVector();           
      // HYPRE_ERROR_ARG
      
      // const auto t0 = std::chrono::system_clock::now();
      // result = HYPRE_AMSSetup(solver, A_parcsr, B_par, x_par);
      // result = HYPRE_AMSSolve(solver, A_parcsr, B_par, x_par);
      HYPRE_BoomerAMGSetup(solver, A_parcsr, B_par, x_par);      
      HYPRE_BoomerAMGSolve(solver, A_parcsr, B_par, x_par);

      // char msg[2048];
      // if (result)
      //   HYPRE_DescribeError(result, msg);
      // HYPRE_ClearError(HYPRE_ERROR_CONV);
      
      // const auto t1 = std::chrono::system_clock::now();
      // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";
      HYPRE_BoomerAMGDestroy(solver);                  

      // result = HYPRE_AMSDestroy(solver);    
      
      x_vector.GetValues(indices_count, indices_ptr, output_ptr); 
    }
    
    void Solve(int count, HYPRE_Real* output_ptr) {
      std::vector<HYPRE_Int> x_indices(count);
      for (HYPRE_Int i = 0; i < count; ++i)
        x_indices[i] = i;
      
      Solve(x_indices.data(), count, output_ptr);  
    }

    auto Solve() {
      std::vector<double> x(nrows);
      Solve(nrows, x.data());
      return x;
    }
  };


  
}
