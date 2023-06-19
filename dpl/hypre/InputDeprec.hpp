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


#include <dpl/general.hpp>
#include <dpl/hypre/ij_matrix.hpp>
#include <dpl/hypre/sparse_matrix_builder.hpp>
#include <dpl/soa.hpp>

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
  namespace keys {
    def_static_key(index)
    def_static_key(value)
  }

  inline constexpr struct populate_tag {} populate;

  /**
   * \brief
   *    nrows = length(ncols) = length(b), number of rows
   *    sum(ncols) = length(cols) = length(values), number of non-zero coefficients
   */
  struct ls_known_ref
  {
    HYPRE_BigInt nrows;     // number of rows   
    
    HYPRE_Int* ncols;       // number of non-zero columns
    HYPRE_BigInt* cols;     // non-zero columns
    HYPRE_Complex* values;  // non-zero coefficients

    HYPRE_Complex* b;       // constant terms
  };


#ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
  inline void load(const boost::interprocess::mapped_region& region, ls_known_ref& lkr, std::pair<HYPRE_BigInt, HYPRE_BigInt>*& rows) {
    using namespace boost::interprocess;

    // mapped_region region(smo, read_only);

    auto* ptr = region.get_address();

    lkr.nrows = *(HYPRE_BigInt*)ptr;
    ptr = (char*)ptr + sizeof(HYPRE_BigInt);

    auto coefs_count = *(HYPRE_BigInt*)ptr;
    ptr = (char*)ptr + sizeof(HYPRE_BigInt);


    lkr.ncols = (HYPRE_BigInt*)ptr;
    ptr = (char*)ptr + lkr.nrows * sizeof(HYPRE_BigInt);

    lkr.b = (HYPRE_Complex*)ptr;
    ptr = (char*)ptr + lkr.nrows * sizeof(HYPRE_Complex);

    lkr.cols = (HYPRE_BigInt*)ptr;
    ptr = (char*)ptr + coefs_count * sizeof(HYPRE_BigInt);

    lkr.values = (HYPRE_Complex*)ptr;
    ptr = (char*)ptr + coefs_count * sizeof(HYPRE_Complex);


    #ifdef HYPRE_SEQUENTIAL
      rows.first = 0;
      rows.second = lkr.nrows - 1;
    #else
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      rows = static_cast<std::pair<HYPRE_BigInt, HYPRE_BigInt>*>(ptr) + rank;
    #endif
  }
#endif





  struct ls_unknown_ref
  {
    HYPRE_Int nvalues;
    HYPRE_BigInt* indices;
    HYPRE_Complex* values;
  };

  struct ls_unknown_storage
  {
    dpl::soa<
      keys::index_t, HYPRE_BigInt,
      keys::value_t, HYPRE_Complex
    > data;

    // ls_unknown_storage(HYPRE_BigInt size) {
    //   data.resize(size);
    // }



    // ls_unknown_storage(HYPRE_BigInt nrows, populate_tag) {
    //   #ifdef HYPRE_SEQUENTIAL
    //     static constexpr auto jlower = 0;
    //     const auto jupper = nrows - 1;
    //   #else
    //     auto [jlower, jupper] = mpi_part(nrows);
    //   #endif
    //
    //   auto count = jupper - jlower + 1;
    //
    //   // ls_unknown_storage lus{count}; 
    //
    //   data.resize(count);
    //   for (HYPRE_BigInt i = 0; i < count; ++i)
    //     data[keys::index][i] = jlower + i;
    //
    //   // #ifdef HYPRE_SEQUENTIAL
    //   //   static constexpr auto jlower = 0;
    //   //   const auto jupper = nrows - 1;
    //   // #else
    //   //   auto [jlower, jupper] = mpi_part(size);
    //   // #endif
    //   //
    //   // data.resize(size);
    //   // for (HYPRE_BigInt i = 0; i < size; ++i)
    //   //   data[keys::index][i] = jlower + i;
    // }

    auto get_ref() {
      return ls_unknown_ref {
        static_cast<HYPRE_Int>(data.size()),
        data[keys::index].get(),
        data[keys::value].get()
      };
    }
  };

  



  inline void solve(const ls_known_ref& in, const ls_unknown_ref& out) {                  
      HYPRE_Solver solver;
      
      HYPRE_BoomerAMGCreate(&solver);      

      HYPRE_BoomerAMGSetTol(solver, 1.e-20);
      // HYPRE_BoomerAMGSetCoarsenType(solver, 10);
      // HYPRE_BoomerAMGSetRestriction(solver, 2);
      // HYPRE_BoomerAMGSetTruncFactor(solver, 4);
      // HYPRE_BoomerAMGSetInterpType(solver, 6);
      // HYPRE_BoomerAMGSetStrongThreshold(solver, 0);
      // HYPRE_BoomerAMGSetTruncFactor(solver, 0);

      // ReSharper disable CppInconsistentNaming
      ij_matrix A_matrix{in.nrows, in.ncols, in.cols, in.values};
      ij_vector b_vector{in.nrows, in.b};
      ij_vector x_vector{in.nrows};
      
      auto* ref_A = A_matrix.par_ref();
      auto* ref_b = b_vector.par_ref();
      auto* ref_x = x_vector.par_ref();
      // ReSharper restore CppInconsistentNaming
      
      // const auto t0 = std::chrono::system_clock::now();

      HYPRE_BoomerAMGSetup(solver, ref_A, ref_b, ref_x);
      HYPRE_BoomerAMGSolve(solver, ref_A, ref_b, ref_x);

      // const auto t1 = std::chrono::system_clock::now();
      // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";
      HYPRE_BoomerAMGDestroy(solver);                  
    
      x_vector.get_values(out.nvalues, out.indices, out.values);
    }




  
  struct InputDeprec
  {
    HYPRE_BigInt nrows;    
    std::vector<HYPRE_Int> ncols_per_row; // size of nrows       
    std::vector<HYPRE_Complex> constants; // size of nrows

    std::vector<HYPRE_BigInt> cols_of_coefs; 
    std::vector<HYPRE_Complex> coefs;

    auto get_ref() {
      return ls_known_ref{
        nrows,
        ncols_per_row.data(),
        cols_of_coefs.data(),
        coefs.data(),
        constants.data()
      };
    }

    InputDeprec() = default;

    InputDeprec(const InputDeprec& other) = default;
    InputDeprec(InputDeprec&& other) noexcept = default;
    InputDeprec& operator=(const InputDeprec& other) = default;
    InputDeprec& operator=(InputDeprec&& other) noexcept = default;
    

    InputDeprec(sparse_matrix_builder& m, std::vector<double>&& b) {             
      nrows = m.nrows;
      ncols_per_row.assign(m.nrows, 1);
      constants = std::move(b);            
      
      const auto coef_count = nrows + m.off_diag_count_;
      
      cols_of_coefs.resize(coef_count);
      coefs.assign(coef_count, 0);

            
      for (HYPRE_Int row_idx = 0, coef_idx = 0; row_idx < nrows; row_idx++) {                
        auto diag_idx = coef_idx++;
        cols_of_coefs[diag_idx] = row_idx;
        coefs[diag_idx] = m.diag[row_idx];        
        
        for (auto [j, j_off_coef] : m.off_diag[row_idx]) {
          ++ncols_per_row[row_idx];
          
          auto off_diag_idx = coef_idx++;
          cols_of_coefs[off_diag_idx] = j;          
          coefs[off_diag_idx] = j_off_coef;          
        }                                                
      }          
    }
    
    // void Solve(HYPRE_BigInt* indices_ptr, HYPRE_Int indices_count, HYPRE_Complex* output_ptr) {                  
    //   HYPRE_Solver solver;
    //   
    //   HYPRE_BoomerAMGCreate(&solver);      
    //
    //   HYPRE_BoomerAMGSetTol(solver, 1.e-20);
    //   
    //   // ReSharper disable once CppInconsistentNaming
    //   ij_matrix A_ij_matrix{nrows, ncols_per_row.data(), cols_of_coefs.data(), coefs.data()};
    //   ij_vector b_ij_vector{nrows, constants.data()};
    //   ij_vector x_ij_vector{nrows};
    //   
    //   auto* ref_A = A_ij_matrix.par_ref();
    //   auto* ref_b = b_ij_vector.par_ref();
    //   auto* ref_x = x_ij_vector.par_ref();           
    //   
    //   // const auto t0 = std::chrono::system_clock::now();
    //   // result = HYPRE_AMSSetup(solver, A_parcsr, B_par, x_par);
    //   // result = HYPRE_AMSSolve(solver, A_parcsr, B_par, x_par);
    //
    //   HYPRE_BoomerAMGSetup(solver, ref_A, ref_b, ref_x);
    //   HYPRE_BoomerAMGSolve(solver, ref_A, ref_b, ref_x);
    //
    //   // char msg[2048];
    //   // if (result)
    //   //   HYPRE_DescribeError(result, msg);
    //   // HYPRE_ClearError(HYPRE_ERROR_CONV);
    //   
    //   // const auto t1 = std::chrono::system_clock::now();
    //   // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";
    //   HYPRE_BoomerAMGDestroy(solver);                  
    //
    //   x_ij_vector.get_values(indices_count, indices_ptr, output_ptr);
    // }
    
    // void Solve(HYPRE_BigInt count, HYPRE_Complex* output_ptr) {
    //   std::vector<HYPRE_Int> indices(count);
    //   for (HYPRE_Int i = 0; i < count; ++i)
    //     indices[i] = i;
    //   
    //   Solve(indices.data(), count, output_ptr);  
    // }

    // auto Solve() {
    //   #ifdef HYPRE_SEQUENTIAL
    //     static constexpr auto jlower = 0;
    //     const auto jupper = nrows - 1;
    //   #else
    //     auto [jlower, jupper] = mpi_part(nrows);
    //   #endif
    //
    //   auto count = jupper - jlower + 1;
    //
    //   ls_unknown_storage lus{count, jlower}; 
    //
    //   // for (HYPRE_BigInt i = 0; i < count; ++i)
    //   //   lus.data[keys::index][i] = jlower + i;
    //
    //   solve(get_ref(), lus.get_ref());
    //   
    //   std::vector<double> x(count);
    //   for (HYPRE_BigInt i = 0; i < count; ++i)
    //     x[i] = lus.data[keys::value][i];
    //
    //   return x;
    // }

  };




  #ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
    using smo_t = boost::interprocess::shared_memory_object;

    inline void save(const InputDeprec& input, const std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>>& blocks, smo_t& smo) {
      // boost::interprocess::shared_memory_object shm (create_only, "MySharedMemory", read_write);
      using namespace boost::interprocess;
      
      auto coefs_count = input.cols_of_coefs.size();
      
      auto buffer_size = sizeof(HYPRE_BigInt) + sizeof(HYPRE_BigInt) +
        input.nrows*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex))
      + coefs_count*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
        blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>);
      
      smo.truncate(buffer_size);
      
      mapped_region region(smo, read_write);
      auto* ptr = region.get_address();

      *(HYPRE_BigInt*)ptr = input.nrows;
      ptr = (char*)ptr + sizeof(HYPRE_BigInt);
      
      *(HYPRE_BigInt*)ptr = coefs_count;
      ptr = (char*)ptr + sizeof(HYPRE_BigInt);
      
      std::memcpy(ptr, input.ncols_per_row.data(), input.nrows*sizeof(HYPRE_BigInt));
      ptr = (char*)ptr + input.nrows*sizeof(HYPRE_BigInt);

      std::memcpy(ptr, input.constants.data(), input.nrows*sizeof(HYPRE_Complex));
      ptr = (char*)ptr + input.nrows*sizeof(HYPRE_Complex);

      std::memcpy(ptr, input.cols_of_coefs.data(), coefs_count*sizeof(HYPRE_BigInt));
      ptr = (char*)ptr + coefs_count*sizeof(HYPRE_BigInt);

      std::memcpy(ptr, input.coefs.data(), coefs_count*sizeof(HYPRE_Complex));
      ptr = (char*)ptr + coefs_count*sizeof(HYPRE_Complex);

      std::memcpy(ptr, blocks.data(), blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>));
    }

    // inline void load(smo_t& smo, InputDeprec& input) {
    //   using namespace boost::interprocess;
    //
    //   mapped_region region(smo, read_only);
    //
    //   auto* ptr = region.get_address();
    //
    //   input.nrows = *(HYPRE_BigInt*)ptr;
    //   ptr = (char*)ptr + sizeof(HYPRE_BigInt);
    //
    //   auto coefs_count = *(HYPRE_BigInt*)ptr;
    //   ptr = (char*)ptr + sizeof(HYPRE_BigInt);
    //
    //
    //   input.ncols_per_row.resize(input.nrows);
    //   std::memcpy(input.ncols_per_row.data(), ptr, input.nrows*sizeof(HYPRE_BigInt));
    //   ptr = (char*)ptr + input.nrows*sizeof(HYPRE_BigInt);
    //   
    //   input.constants.resize(input.nrows);
    //   std::memcpy(input.constants.data(), ptr, input.nrows*sizeof(HYPRE_Complex));
    //   ptr = (char*)ptr + input.nrows*sizeof(HYPRE_Complex);
    //   
    //   input.cols_of_coefs.resize(coefs_count);
    //   std::memcpy(input.cols_of_coefs.data(), ptr, coefs_count*sizeof(HYPRE_BigInt));
    //   ptr = (char*)ptr + coefs_count*sizeof(HYPRE_BigInt);
    //   
    //   input.coefs.resize(coefs_count);
    //   std::memcpy(input.coefs.data(), ptr, coefs_count*sizeof(HYPRE_Complex));
    //   // ptr += nrows*sizeof(HYPRE_Complex);
    // }
#endif

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