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

namespace dpl::hypre
{
  /**
   * \brief
   *    nrows = length(ncols) = length(b), number of rows
   *    nvalues = sum(ncols) = length(cols) = length(values), number of non-zero coefficients
   */
  struct ls_known_ref
  {
    HYPRE_Int* ncols;       // number of non-zero columns
    HYPRE_BigInt* cols;     // non-zero columns
    HYPRE_Complex* values;  // non-zero coefs

    HYPRE_Complex* b;       // constant terms
  };

  struct ls_known_storage
  {
    std::unique_ptr<HYPRE_Int[]> ncols;       
    std::unique_ptr<HYPRE_BigInt[]> cols;     
    std::unique_ptr<HYPRE_Complex[]> values;  

    std::unique_ptr<HYPRE_Complex[]> b;       

    operator ls_known_ref() const {
      return { ncols.get(), cols.get(), values.get(), b.get() };
    }
  };

  class ls_known_storage_builder
  {
    HYPRE_BigInt nrows_;
    size_t nvalues_;
    ls_known_storage lks_;

    std::unique_ptr<size_t[]> diag_shift_;
    std::unique_ptr<HYPRE_Int[]> off_relative_;
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



  inline std::pair<HYPRE_Real, HYPRE_Int> solve(
    HYPRE_BigInt ilower, HYPRE_BigInt iupper, const ls_known_ref& in, HYPRE_Complex* values,
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

    auto nrows = iupper - ilower + 1;
    auto indices = std::make_unique<HYPRE_BigInt[]>(nrows);
    std::iota(indices.get(), indices.get() + nrows, ilower);

    ij_matrix A{ilower, iupper, indices.get(), in.ncols, in.cols, in.values};
    ij_vector b{ilower, iupper, in.b};
    ij_vector x{ilower, iupper};



    // const auto t0 = std::chrono::system_clock::now();

    // int w_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    //
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_0\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    HYPRE_BoomerAMGSetup(solver, A, b, x);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_1\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);
    
    HYPRE_BoomerAMGSolve(solver, A, b, x);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == 0)
    //   std::cout << "\nPOS_2\n" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);


    // const auto t1 = std::chrono::system_clock::now();
    // std::cout << "Actual setup&solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms\n";


    HYPRE_Real final_residual = 0;
    HYPRE_Int num_interations = 0;

    #ifndef HYPRE_SEQUENTIAL
    int w_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    if (w_rank == 0)
    #endif
    {
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_residual);
      HYPRE_BoomerAMGGetNumIterations(solver, &num_interations);
    }

    HYPRE_BoomerAMGDestroy(solver);                  
  
    x.get_values(nrows, indices.get(), values);

    return {final_residual, num_interations};
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
}


#ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace dpl::hypre
{
  namespace mpi
  {
    inline constexpr auto smo_name = "xpm-hypre-input";
    namespace bi = boost::interprocess;

    class block_info
    {
      bi::mapped_region region_;
      ls_known_ref lkr_;

    public:
      explicit block_info(int rank)
        : region_{bi::shared_memory_object{bi::open_only, smo_name, bi::read_only}, bi::read_only} {

        size_t nvalues;  
        
        parser p{region_.get_address()};
        
        p.read(global_nrows);
        p.read(nvalues);
        p.read(lkr_.ncols, global_nrows);
        p.read(lkr_.b, global_nrows);
        p.read(lkr_.cols, nvalues);
        p.read(lkr_.values, nvalues);

        local_rows = static_cast<std::pair<HYPRE_BigInt, HYPRE_BigInt>*>(p.ptr()) + rank;
      }

      HYPRE_BigInt global_nrows;
      std::pair<HYPRE_BigInt, HYPRE_BigInt>* local_rows;

      operator ls_known_ref() const { return lkr_; }
    };

    inline auto load(int rank) {
      return block_info{rank};
    }

    inline void save(
      const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>>& blocks) {

      bi::shared_memory_object smo{bi::open_or_create, smo_name, bi::read_write};

      auto buffer_size = sizeof(HYPRE_BigInt) + sizeof(size_t) +
        nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
        nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
        blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>);
      
      smo.truncate(buffer_size);  // NOLINT(cppcoreguidelines-narrowing-conversions)
      
      bi::mapped_region region(smo, bi::read_write);

      parser p{region.get_address()};
      p.write(nrows);
      p.write(nvalues);
      p.write(input.ncols, nrows);
      p.write(input.b, nrows);
      p.write(input.cols, nvalues);
      p.write(input.values, nvalues);
      
      std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(std::pair<HYPRE_BigInt, HYPRE_BigInt>));
    }
  }
}
#endif