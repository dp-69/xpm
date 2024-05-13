/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "core.hpp"

#include <vector>

namespace dpl::hypre
{
  class sparse_matrix_builder
  {                    
    HYPRE_BigInt off_diag_count_;    
    
    friend struct InputDeprec;

  public:
    HYPRE_BigInt nrows;
    
    std::vector<HYPRE_Complex> diag;
    std::vector<std::forward_list<std::tuple<HYPRE_BigInt, HYPRE_Complex>>> off_diag;
    
    sparse_matrix_builder() = default;
    
    explicit sparse_matrix_builder(HYPRE_BigInt n) {
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
    

    InputDeprec(sparse_matrix_builder& m, std::vector<HYPRE_Complex>&& b) {             
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


// #ifdef DPL_HYPRE_BOOST_SHARED_MEMORY
// #include <boost/interprocess/shared_memory_object.hpp>
// #include <boost/interprocess/mapped_region.hpp>
//
// namespace dpl::hypre
// {
//   namespace mpi
//   {
//     class parser
//     {
//       void* ptr_;
//
//     public:
//       explicit parser(void* ptr) : ptr_(ptr) {}
//
//       template <typename T>
//       void read(T& val) {
//         val = *static_cast<T*>(ptr_);
//         ptr_ = static_cast<char*>(ptr_) + sizeof(T); 
//       }
//       
//       template <typename T, typename S>
//       void read(T*& val, S size) {
//         val = static_cast<T*>(ptr_);
//         ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
//       }
//
//       template <typename T>
//       void write(T val) {
//         *static_cast<T*>(ptr_) = val;
//         ptr_ = static_cast<char*>(ptr_) + sizeof(T);
//       }
//
//       template <typename T, typename S>
//       void write(T* val, S size) {
//         std::memcpy(ptr_, val, size*sizeof(T));
//         ptr_ = static_cast<char*>(ptr_) + size*sizeof(T);
//       }
//
//       auto* ptr() const { return ptr_; }
//     };
//
//     inline constexpr auto smo_hypre_input = "hypre-input";
//     inline constexpr auto smo_hypre_output = "hypre-output";
//
//     namespace bi = boost::interprocess;
//
//     class block_info
//     {
//       bi::mapped_region region_;
//       ls_known_ref lkr_;
//
//     public:
//       explicit block_info(int rank)
//         : region_{bi::shared_memory_object{bi::open_only, smo_hypre_input, bi::read_only}, bi::read_only} {
//
//         size_t nvalues;  
//         
//         parser p{region_.get_address()};
//         
//         p.read(global_nrows);
//         p.read(nvalues);
//         p.read(lkr_.ncols, global_nrows);
//         p.read(lkr_.b, global_nrows);
//         p.read(lkr_.cols, nvalues);
//         p.read(lkr_.values, nvalues);
//
//         range = static_cast<index_range*>(p.ptr()) + rank;
//       }
//
//       HYPRE_BigInt global_nrows;
//       index_range* range;
//
//       operator ls_known_ref() const { return lkr_; }
//     };
//
//     inline auto load_values(HYPRE_BigInt nrows) {
//       auto values = std::make_unique<HYPRE_Complex[]>(nrows);
//       std::memcpy(
//         values.get(),
//         bi::mapped_region{
//           bi::shared_memory_object{bi::open_only, smo_hypre_output, bi::read_only},
//           bi::read_only}.get_address(),
//         nrows*sizeof(HYPRE_Complex));
//       return values;
//     }
//
//     inline void save_values(const HYPRE_Complex* values, HYPRE_BigInt nrows) {
//       bi::shared_memory_object smo_output{bi::open_or_create, smo_hypre_output, bi::read_write};
//       smo_output.truncate(nrows*sizeof(HYPRE_Complex));  // NOLINT(cppcoreguidelines-narrowing-conversions)
//       bi::mapped_region region_output(smo_output, bi::read_write);
//       std::memcpy(region_output.get_address(), values, nrows*sizeof(HYPRE_Complex));
//     }
//
//     inline auto load_block(int rank) { return block_info{rank}; }
//
//     inline void save(
//       const ls_known_ref& input, HYPRE_BigInt nrows, size_t nvalues, const std::vector<index_range>& blocks) {
//
//       bi::shared_memory_object smo{bi::open_or_create, smo_hypre_input, bi::read_write};
//
//       auto buffer_size = sizeof(HYPRE_BigInt) + sizeof(size_t) +
//         nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
//         nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)) +
//         blocks.size()*sizeof(index_range);
//       
//       smo.truncate(buffer_size);  // NOLINT(cppcoreguidelines-narrowing-conversions)
//       
//       bi::mapped_region region(smo, bi::read_write);
//
//       parser p{region.get_address()};
//       p.write(nrows);
//       p.write(nvalues);
//       p.write(input.ncols, nrows);
//       p.write(input.b, nrows);
//       p.write(input.cols, nvalues);
//       p.write(input.values, nvalues);
//       
//       std::memcpy(p.ptr(), blocks.data(), blocks.size()*sizeof(index_range));
//     }
//   }
// }
// #endif