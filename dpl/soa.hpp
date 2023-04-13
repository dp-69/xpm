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
#include "dpl/static_map.hpp"

#include <memory>
#include <vector>


namespace dpl
{
  // template<typename T, typename = void>
  // struct static_item_desc {
  //   using key_type = T;
  //   using value_type = std::vector<typename T::type>;
  // };
  //
  // template <typename T>
  // struct static_item_desc<T, std::void_t<typename T::value_type>> {
  //   using key_type = T;
  //   using value_type = typename T::value_type;
  // };
  /*static_item_desc<Ts>...*/

  // template <typename... Args>
  // class attributes_soa_uptr_impl_;
  //
  // template <>
  // class attributes_soa_uptr_impl_<>
  // {
  // };
  //
  // template <typename... Args>
  // class attributes_soa_uptr : public static_map<Args...>
  // {
  //   
  // };
  

  template <typename... Args>
  class soa_impl_ : public static_map_impl_<Args...>
  {
    using size_type = std::size_t;

    template <typename T>
    static void resize_(std::unique_ptr<T[]>& uptr, size_type size) {
      uptr = std::make_unique<T[]>(size);
    }

    template <typename T>
    static void assign_(std::unique_ptr<T[]>& uptr, size_type size, const T& v) {
      resize_(uptr, size);
      for (size_type i = 0; i < size; ++i)
        uptr[i] = v;
    }

    // template <typename T>
    // static void resize_(std::vector<T>& vec, size_type size) {
    //   vec.resize(size);
    // }
    //
    // template <typename T>
    // static void assign_(std::vector<T>& vec, size_type size, const T& v) {
    //   vec.assign(size, v);
    // }
    
  public:
    soa_impl_() = default;

    explicit soa_impl_(Args...) {}
    
    void resize(const size_type count) {
      this->apply([=](auto& attrib) { resize_(attrib, count); });
    }

    template <typename T>
    void assign(const size_type count, const T& v) {
      this->apply([=](auto& attrib) { assign_(attrib, count, v); });
    }
  };


  namespace type_map
  {
    struct unique_ptr_array
    {
      template <typename T>
      using type = std::unique_ptr<T[]>;
    };
  }
  
  template <typename... Args>
  class soa : public soa_impl_<type_map::unique_ptr_array, Args...> {};
}