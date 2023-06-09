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

#include <dpl/hypre/ij_matrix.hpp>
#include <dpl/hypre/sparse_matrix.hpp>

#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif

// #include "rrm/fd/core/linear_solver/HypreVectorMatrix.hpp"

// #include <boost/format/format_fwd.hpp>
// #include <boost/format.hpp>

#include <vector>
#include <fstream>

#include <forward_list>
#include <iostream>

#ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#endif


namespace dpl::hypre
{
  



  struct InputPtr
  {
    HYPRE_Int nrows;    
    HYPRE_Int* ncols; // size of nrows
    HYPRE_Complex* b; // constants, b vector, size of nrows

    HYPRE_BigInt* cols; 
    HYPRE_Complex* values;
  };


  
  struct Input
  {
    HYPRE_Int nrows;    
    std::vector<HYPRE_Int> ncols_per_row; // size of nrows       
    std::vector<HYPRE_Real> constants; // size of nrows

    std::vector<HYPRE_Int> cols_of_coefs; 
    std::vector<HYPRE_Real> coefs;

#ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
    using smo_t = boost::interprocess::shared_memory_object;
    
    void Save(smo_t& smo) {
      // boost::interprocess::shared_memory_object shm (create_only, "MySharedMemory", read_write);
      using namespace boost::interprocess;
      
      auto coefs_count = cols_of_coefs.size();
      
      auto buffer_size = sizeof(HYPRE_Int) + sizeof(HYPRE_Int) +
        nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Real))
      + coefs_count*(sizeof(HYPRE_Int) + sizeof(HYPRE_Real));
      
      smo.truncate(buffer_size);
      
      mapped_region region(smo, read_write);
      auto* ptr = region.get_address();

      *(HYPRE_Int*)ptr = nrows;
      ptr = (char*)ptr + sizeof(HYPRE_Int);
      
      *(HYPRE_Int*)ptr = coefs_count;
      ptr = (char*)ptr + sizeof(HYPRE_Int);
      
      std::memcpy(ptr, ncols_per_row.data(), nrows*sizeof(HYPRE_Int));
      ptr = (char*)ptr + nrows*sizeof(HYPRE_Int);

      std::memcpy(ptr, constants.data(), nrows*sizeof(HYPRE_Real));
      ptr = (char*)ptr + nrows*sizeof(HYPRE_Real);

      std::memcpy(ptr, cols_of_coefs.data(), coefs_count*sizeof(HYPRE_Int));
      ptr = (char*)ptr + coefs_count*sizeof(HYPRE_Int);

      std::memcpy(ptr, coefs.data(), coefs_count*sizeof(HYPRE_Real));
      // ptr += nrows*sizeof(HYPRE_Real);
    }

    void Load(smo_t& smo) {
      using namespace boost::interprocess;

      mapped_region region(smo, read_only);

      auto* ptr = region.get_address();

      nrows = *(HYPRE_Int*)ptr;
      ptr = (char*)ptr + sizeof(HYPRE_Int);

      auto coefs_count = *(HYPRE_Int*)ptr;
      ptr = (char*)ptr + sizeof(HYPRE_Int);


      ncols_per_row.resize(nrows);
      std::memcpy(ncols_per_row.data(), ptr, nrows*sizeof(HYPRE_Int));
      ptr = (char*)ptr + nrows*sizeof(HYPRE_Int);
      
      constants.resize(nrows);
      std::memcpy(constants.data(), ptr, nrows*sizeof(HYPRE_Real));
      ptr = (char*)ptr + nrows*sizeof(HYPRE_Real);
      
      cols_of_coefs.resize(coefs_count);
      std::memcpy(cols_of_coefs.data(), ptr, coefs_count*sizeof(HYPRE_Int));
      ptr = (char*)ptr + coefs_count*sizeof(HYPRE_Int);
      
      coefs.resize(coefs_count);
      std::memcpy(coefs.data(), ptr, coefs_count*sizeof(HYPRE_Real));
      // ptr += nrows*sizeof(HYPRE_Real);
    }

#endif

    

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
    
    // void Export(const std::string& filename, bool numerical = false) {
    //   std::ofstream out(filename);      
    //
    //   auto* colPtr = cols_of_coefs.data();
    //   auto* coefPtr = coefs.data();
    //   for (auto i = 0; i < nrows; ++i) {
    //     std::vector<HYPRE_Real> row_values(nrows, 0);
    //
    //     for (auto j = 0; j < ncols_per_row[i]; ++j)
    //       row_values[*colPtr++] = *coefPtr++;
    //
    //     for (auto j = 0; j < nrows; ++j) {          
    //       if (numerical)
    //         out << boost::format("%10.2e") % row_values[j];
    //       else
    //         out << (row_values[j] ? row_values[j] == 1 ? '1' : 'v' : '0'); // NOLINT(clang-diagnostic-float-conversion, clang-diagnostic-float-equal)
    //       
    //       out << ' ';
    //     }
    //
    //     out << "| ";
    //             
    //     if (numerical)
    //       out << boost::format("%10.2e") % constants[i];
    //     else
    //       out << (constants[i] ? constants[i] == 1 ? '1' : 'v' : '0');  // NOLINT(clang-diagnostic-float-conversion, clang-diagnostic-float-equal)
    //     
    //     out << std::endl;
    //   }
    //   
    //   // out << "one" << " two" << boost::format("%e %.3e") % 1.5e3 % 0.0048799;
    //
    //   out.close();
    // }
    
    void Solve(HYPRE_Int* indices_ptr, HYPRE_Int indices_count, HYPRE_Real* output_ptr) {                  
      HYPRE_Solver solver;
      
      HYPRE_BoomerAMGCreate(&solver);      

      HYPRE_BoomerAMGSetTol(solver, 1.e-20);
      
      // ReSharper disable once CppInconsistentNaming
      ij_matrix A_ij_matrix{nrows, ncols_per_row.data(), cols_of_coefs.data(), coefs.data()};
      ij_vector b_ij_vector{nrows, constants.data()};
      ij_vector x_ij_vector{nrows};
      
      auto* ref_A = A_ij_matrix.par_ref();
      auto* ref_b = b_ij_vector.par_ref();
      auto* ref_x = x_ij_vector.par_ref();           
      
      // const auto t0 = std::chrono::system_clock::now();
      // result = HYPRE_AMSSetup(solver, A_parcsr, B_par, x_par);
      // result = HYPRE_AMSSolve(solver, A_parcsr, B_par, x_par);

      HYPRE_BoomerAMGSetup(solver, ref_A, ref_b, ref_x);
      HYPRE_BoomerAMGSolve(solver, ref_A, ref_b, ref_x);

      // char msg[2048];
      // if (result)
      //   HYPRE_DescribeError(result, msg);
      // HYPRE_ClearError(HYPRE_ERROR_CONV);
      
      // const auto t1 = std::chrono::system_clock::now();
      // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";
      HYPRE_BoomerAMGDestroy(solver);                  

      x_ij_vector.get_values(indices_count, indices_ptr, output_ptr);
    }
    
    void Solve(HYPRE_Int count, HYPRE_Real* output_ptr) {
      std::vector<HYPRE_Int> indices(count);
      for (HYPRE_Int i = 0; i < count; ++i)
        indices[i] = i;
      
      Solve(indices.data(), count, output_ptr);  
    }

    auto Solve() {
#ifndef HYPRE_SEQUENTIAL
      auto [ilower, iupper] = mpi_part(nrows);

      auto count = iupper - ilower + 1;

      std::vector<double> x(count);
      std::vector<HYPRE_Int> indices(count);
      for (HYPRE_BigInt i = 0; i < count; ++i)
        indices[i] = ilower + i;
      
      Solve(indices.data(), count, x.data());
      
      return x;
#else
      std::vector<double> x(nrows);
      Solve(nrows, x.data());
      return x;
#endif
    }
  };


  
}




 // (Optional) Set the convergence tolerance, if BoomerAMG is used as a solver. If it is used as a preconditioner,
      //   it should be set to 0. The default is 1.e-7.

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