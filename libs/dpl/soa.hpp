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
#include <dpl/static_map.hpp>

#include <memory>
#include <vector>

#ifdef __cpp_lib_span
  #include <span>    
#endif


namespace dpl
{
  // template <typename T, typename = void>
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

    size_type size_ = 0;
    
    template <typename T>
    void resize_(std::unique_ptr<T[]>& uptr, size_type new_size) {
      auto new_uptr = std::make_unique<T[]>(new_size);
      std::copy(uptr.get(), uptr.get() + std::min(size_, new_size), new_uptr.get());
      uptr = std::move(new_uptr);
    }

    template <typename T>
    void assign_(std::unique_ptr<T[]>& uptr, size_type size, const T& v) {
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

    auto size() const {
      return size_;
    }

    void resize(const size_type count) {
      this->apply([this, count](auto& attrib) {
        this->resize_(attrib, count);
      });
      size_ = count;
    }

    template <typename T>
    void assign(const size_type count, const T& v) {
      this->apply([=](auto& attrib) { assign_(attrib, count, v); }); // TODO, tidy up
      size_ = count;
    }

    template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    const auto* ptr(Key k) const {
      static constexpr auto lookup = static_map_rec_<Args...>::find(k);
      static_assert(lookup, "Type key is not available");
      return lookup(this).get();
    }

    template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    auto* ptr(Key k) {
      static constexpr auto lookup = static_map_rec_<Args...>::find(k);
      static_assert(lookup, "Type key is not available");
      return lookup(this).get();
    }

    #ifdef __cpp_lib_span
    template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    auto span(Key k) const {
      static constexpr auto lookup = static_map_rec_<Args...>::find(k);
      static_assert(lookup, "Type key is not available");
      auto* ptr = lookup(this).get();
      return std::span{ptr, ptr + size_};
      // return std::ranges::subrange{ptr, ptr + size_};
    }    
    #endif
    
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

  // template <typename... Args>
  // class soa_size : public soa_impl_<type_map::unique_ptr_array, Args...> {};
}