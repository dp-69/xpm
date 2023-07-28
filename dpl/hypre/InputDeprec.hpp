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


#include <dpl/general.hpp>
#include <dpl/hypre/ij_matrix.hpp>
#include <dpl/hypre/sparse_matrix_builder.hpp>
#include <dpl/soa.hpp>

#include <vector>

#ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#endif


namespace dpl::hypre
{
  // namespace keys {
  //   def_static_key(index)
  //   def_static_key(value)
  // }

  // inline constexpr struct populate_tag {} populate;

  /**
   * \brief
   *    nrows = length(ncols) = length(b), number of rows
   *    nvalues = sum(ncols) = length(cols) = length(values), number of non-zero coefficients
   */
  struct ls_known_ref
  {
    // HYPRE_BigInt nrows;     // number of rows   
    
    HYPRE_Int* ncols;       // number of non-zero columns
    HYPRE_BigInt* cols;     // non-zero columns
    HYPRE_Complex* values;  // non-zero coefs

    HYPRE_Complex* b;       // constant terms
  };

  struct ls_known_storage
  {
    // HYPRE_BigInt nrows;                       // number of rows

    std::unique_ptr<HYPRE_Int[]> ncols;       // number of non-zero columns
    std::unique_ptr<HYPRE_Complex[]> b;       // constant terms

    // size_t nvalues;                           // number of non-zero coefs
    std::unique_ptr<HYPRE_BigInt[]> cols;     // non-zero columns
    std::unique_ptr<HYPRE_Complex[]> values;  // non-zero coefs

    auto get_ref() {
      return ls_known_ref {
        // nrows,
        ncols.get(),
        cols.get(),
        values.get(),
        b.get()
      };
    }
  };

  class ls_known_storage_builder
  {
    ls_known_storage lks_;

    std::unique_ptr<size_t[]> diag_shift_;
    std::unique_ptr<HYPRE_Int[]> off_relative_;

    HYPRE_BigInt nrows_;
    size_t nvalues_;

  public:
    auto nvalues() const { return nvalues_; }

    void allocate_rows(HYPRE_BigInt nrows) {
      nrows_ = nrows;
      lks_.ncols = std::make_unique<HYPRE_Int[]>(nrows);
      lks_.b = std::make_unique<HYPRE_Complex[]>(nrows);

      for (HYPRE_BigInt i = 0; i < nrows; ++i) {
        lks_.ncols[i] = 1; // diag coef
        lks_.b[i] = 0;
      }

      nvalues_ = nrows;
    }

    void reserve_connection(HYPRE_BigInt i0, HYPRE_BigInt i1) {
      ++lks_.ncols[i0];
      ++lks_.ncols[i1];
      nvalues_ += 2;
    }

    void allocate_values() {
      lks_.cols = std::make_unique<HYPRE_BigInt[]>(nvalues_);
      lks_.values = std::make_unique<HYPRE_Complex[]>(nvalues_);

      diag_shift_ = std::make_unique<size_t[]>(nrows_);
      off_relative_ = std::make_unique<HYPRE_Int[]>(nrows_);

      diag_shift_[0] = 0;
      off_relative_[0] = 0; // pre increment will be required
      lks_.cols[0] = 0;
      lks_.values[0] = 0; // diag
      
      for (HYPRE_BigInt i = 1; i < nrows_; ++i) {
        auto s = diag_shift_[i - 1] + lks_.ncols[i - 1];
        diag_shift_[i] = s;
        off_relative_[i] = 0; // pre increment will be required
        lks_.cols[s] = i;
        lks_.values[s] = 0; // diag
      }
    }

    void add_b(HYPRE_BigInt i, double value) {
      lks_.b[i] += value;
    }

    void add_diag(HYPRE_BigInt i, double value) {
      lks_.values[diag_shift_[i]] += value;
    }

    void set_connection(HYPRE_BigInt i0, HYPRE_BigInt i1, double value) {
      auto shift0 = diag_shift_[i0];
      lks_.values[shift0] += value; // diag
      shift0 += ++off_relative_[i0];
      lks_.cols[shift0] = i1;
      lks_.values[shift0] = -value;

      auto shift1 = diag_shift_[i1];
      lks_.values[shift1] += value; // diag
      shift1 += ++off_relative_[i1];
      lks_.cols[shift1] = i0;
      lks_.values[shift1] = -value;
    }

    auto acquire_storage() {
      return std::move(lks_);
    }
  };


  class parser
  {
    void* ptr_;

  public:
    explicit parser(void* ptr) : ptr_(ptr) {}

    template <typename T>
    void read(T& val) {
      val = *static_cast<T*>(ptr_);
      ptr_ = static_cast<char*>(ptr_) + sizeof(T); 
    }
    
    template <typename T, typename S>
    void read(T*& val, S size) {
      val = static_cast<T*>(ptr_);
      ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    }

    // std::memcpy(ptr, input.ncols_per_row.data(), input.nrows*sizeof(HYPRE_Int));
    //   ptr = (char*)ptr + input.nrows*sizeof(HYPRE_Int);

    template <typename T>
    void write(T val) {
      *static_cast<T*>(ptr_) = val;
      ptr_ = static_cast<char*>(ptr_) + sizeof(T);
    }

    template <typename T, typename S>
    void write(T* val, S size) {
      std::memcpy(ptr_, val, size*sizeof(T));
      ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
    }

    auto* ptr() const { return ptr_; }
  };




  // struct ls_unknown_ref
  // {
  //   HYPRE_Int nvalues;
  //   HYPRE_BigInt* indices;
  //   HYPRE_Complex* values;
  // };

  // struct ls_unknown_storage
  // {
  //   dpl::soa<
  //     keys::index_t, HYPRE_BigInt,
  //     keys::value_t, HYPRE_Complex
  //   > data;
  //
  //   explicit ls_unknown_storage(HYPRE_BigInt size, HYPRE_BigInt jlower = 0) {
  //     data.resize(size);
  //     std::iota(data[keys::index].get(), data[keys::index].get() + size, jlower);
  //   }
  //
  //
  //
  //   // ls_unknown_storage(HYPRE_BigInt nrows, populate_tag) {
  //   //   #ifdef HYPRE_SEQUENTIAL
  //   //     static constexpr auto jlower = 0;
  //   //     const auto jupper = nrows - 1;
  //   //   #else
  //   //     auto [jlower, jupper] = mpi_part(nrows);
  //   //   #endif
  //   //
  //   //   auto count = jupper - jlower + 1;
  //   //
  //   //   // ls_unknown_storage lus{count}; 
  //   //
  //   //   data.resize(count);
  //   //   for (HYPRE_BigInt i = 0; i < count; ++i)
  //   //     data[keys::index][i] = jlower + i;
  //   //
  //   //   // #ifdef HYPRE_SEQUENTIAL
  //   //   //   static constexpr auto jlower = 0;
  //   //   //   const auto jupper = nrows - 1;
  //   //   // #else
  //   //   //   auto [jlower, jupper] = mpi_part(size);
  //   //   // #endif
  //   //   //
  //   //   // data.resize(size);
  //   //   // for (HYPRE_BigInt i = 0; i < size; ++i)
  //   //   //   data[keys::index][i] = jlower + i;
  //   // }
  //
  //   auto get_ref() {
  //     return ls_unknown_ref {
  //       static_cast<HYPRE_Int>(data.size()),
  //       data[keys::index].get(),
  //       data[keys::value].get()
  //     };
  //   }
  // };

  



  inline std::pair<HYPRE_Real, HYPRE_Int> solve(HYPRE_BigInt ilower, HYPRE_BigInt iupper, const ls_known_ref& in, HYPRE_Complex* values,
    HYPRE_Real tolerance = 1.e-20, HYPRE_Int max_iterations = 20
    // HYPRE_Real tolerance = 1.e-9, HYPRE_Int max_iterations = 1000
  ) {

    HYPRE_Solver solver;
    
    HYPRE_BoomerAMGCreate(&solver);      

    HYPRE_BoomerAMGSetTol(solver, tolerance);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iterations);

    // HYPRE_BoomerAMGSetMaxLevels(solver, 50);
    // HYPRE_BoomerAMGSetPMaxElmts(solver, 0);
    // HYPRE_BoomerAMGSetMaxCoarseSize(solver, 18);

    // HYPRE_BoomerAMGSetCoarsenType(solver, 0);
    // HYPRE_BoomerAMGSetRestriction(solver, 2);
    // HYPRE_BoomerAMGSetTruncFactor(solver, 4);
    // HYPRE_BoomerAMGSetInterpType(solver, 6);
    // HYPRE_BoomerAMGSetStrongThreshold(solver, 0);
    // HYPRE_BoomerAMGSetTruncFactor(solver, 0);

    

    // #ifdef HYPRE_SEQUENTIAL
    //   static constexpr auto ilower = 0;
    //   const auto iupper = in.nrows - 1;
    // #else
    //   auto [ilower, iupper] = mpi::range; // mpi_part(nrows);
    // #endif

    auto nrows = iupper - ilower + 1;
    auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    std::iota(indices.get(), indices.get() + nrows, ilower);

    // ReSharper disable CppInconsistentNaming
    ij_matrix ij_A{ilower, iupper, indices.get(), in.ncols, in.cols, in.values};
    ij_vector ij_b{ilower, iupper, in.b};
    ij_vector ij_x{ilower, iupper};
    
    auto* ref_A = ij_A.par_ref();
    auto* ref_b = ij_b.par_ref();
    auto* ref_x = ij_x.par_ref();
    // ReSharper restore CppInconsistentNaming




    // if (w_rank == 0) {
    //   HYPRE_BoomerAMGSetLogging(solver, 2);
    // }

    // const auto t0 = std::chrono::system_clock::now();



    // int w_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    //
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_0\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    HYPRE_BoomerAMGSetup(solver, ref_A, ref_b, ref_x);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_1\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);
    
    HYPRE_BoomerAMGSolve(solver, ref_A, ref_b, ref_x);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_2\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);
    //
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);





    // const auto t1 = std::chrono::system_clock::now();
    // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";


    HYPRE_Real final_res = 0;
    HYPRE_Int num_interations = 0;

    {
    #ifndef HYPRE_SEQUENTIAL
    int w_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    if (w_rank == 0)
    #endif
    {
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res);
      HYPRE_BoomerAMGGetNumIterations(solver, &num_interations);
    }
    }
    

    HYPRE_BoomerAMGDestroy(solver);                  
  
    ij_x.get_values(nrows, indices.get(), values);

    return {final_res, num_interations};
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
        // nrows,
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
        cols_of_coefs[diag_idx] = row_idx; // diagonal coef
        coefs[diag_idx] = m.diag[row_idx]; //       
        
        for (auto [j, j_off_coef] : m.off_diag[row_idx]) {
          ++ncols_per_row[row_idx];
          
          auto off_diag_idx = coef_idx++;
          cols_of_coefs[off_diag_idx] = j;          
          coefs[off_diag_idx] = j_off_coef;          
        }                                                
      }          
    }
  };




  #ifdef DPL_HYPRE_BOOST_SHARED_MEMORY

  inline HYPRE_BigInt load(const boost::interprocess::mapped_region& region, ls_known_ref& lkr, std::pair<HYPRE_BigInt, HYPRE_BigInt>*& rows) {
    using namespace boost::interprocess;

    HYPRE_BigInt nrows;
    size_t nvalues;  // NOLINT(cppcoreguidelines-init-variables)
    
    parser p{region.get_address()};
    
    p.read(nrows);
    p.read(nvalues);
    p.read(lkr.ncols, nrows);
    p.read(lkr.b, nrows);
    p.read(lkr.cols, nvalues);
    p.read(lkr.values, nvalues);

    #ifdef HYPRE_SEQUENTIAL
      rows->first = 0;
      rows->second = lkr.nrows - 1;
    #else
      int rank;  // NOLINT(cppcoreguidelines-init-variables)
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      rows = static_cast<std::pair<HYPRE_BigInt, HYPRE_BigInt>*>(p.ptr()) + rank;
    #endif

    return nrows;
  }

  using smo_t = boost::interprocess::shared_memory_object;

  inline void save(const InputDeprec& input, const std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>>& blocks, smo_t& smo) {
    using namespace boost::interprocess;
    
    auto nvalues = input.cols_of_coefs.size();
    
    auto buffer_size = sizeof(HYPRE_BigInt) + sizeof(size_t) +
      input.nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex))
    + nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
      blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>);
    
    smo.truncate(buffer_size);
    
    mapped_region region(smo, read_write);

    parser p{region.get_address()};
    p.write(input.nrows);
    p.write(nvalues);
    p.write(input.ncols_per_row.data(), input.nrows);
    p.write(input.constants.data(), input.nrows);
    p.write(input.cols_of_coefs.data(), nvalues);
    p.write(input.coefs.data(), nvalues);

    std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>));
  }

  inline void save(const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>>& blocks, smo_t& smo) {
    using namespace boost::interprocess;
    
    auto buffer_size = sizeof(HYPRE_BigInt) + sizeof(size_t) +
      nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex))
    + nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
      blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>);
    
    smo.truncate(buffer_size);
    
    mapped_region region(smo, read_write);

    parser p{region.get_address()};
    p.write(nrows);
    p.write(nvalues);
    p.write(input.ncols, nrows);
    p.write(input.b, nrows);
    p.write(input.cols, nvalues);
    p.write(input.values, nvalues);
    
    std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>));
  }
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