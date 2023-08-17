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

#include <boost/pool/object_pool.hpp>

namespace dpl::graph
{
  template <class T,
    typename = std::enable_if_t<sizeof(T) >= sizeof(std::size_t)>>
  class smart_pool
  { 
    static T* get_next(T* x) {       
      return *reinterpret_cast<T**>(x);
    }

    static void set_next(T* x, T* next) {
      *reinterpret_cast<T**>(x) = next;
    }

    boost::object_pool<T> pool_;
    T* released_stack_top_ = nullptr;

  public:
    T* acquire() {
      if (released_stack_top_) {
        auto result = released_stack_top_;        
        released_stack_top_ = get_next(released_stack_top_);       
        return result;
      }      

      return pool_.malloc();
    }

    void release(T* x) {
      set_next(x, released_stack_top_);
      released_stack_top_ = x;
    }
  };
}