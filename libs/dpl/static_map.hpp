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

// #include <dpl/general.hpp>

namespace dpl
{
  namespace type_map
  {
    struct identity
    {
      template <typename T>
      using type = T;
    };
  }
  
  struct lookup_item_tag {};

  template <typename This = void>
  struct lookup_item : lookup_item_tag
  {
    explicit constexpr operator bool() const noexcept {
      return !std::is_same_v<This, void>;
    }

    template <typename T>
    auto& operator()(T* x) const {
      return static_cast<This*>(x)->value_;
    }

    template <typename T>
    const auto& operator()(const T* x) const {
      return static_cast<const This*>(x)->value_;
    }
  };


  template <typename... Args>
  class static_map_rec_;
  
  template <typename ValueMap>
  class static_map_rec_<ValueMap>
  {
  protected:
    template <typename Action>
    static void apply(const Action&) {}

    template <typename Return, typename Pred>
    static Return* first(const Pred&) { return nullptr; }
    
    template <typename Key>
    static constexpr auto find(Key = Key{}) { return lookup_item<>{}; }
  };

  template <typename ValueMap, typename ThisKey, typename ThisValue, typename... Rest>
  class static_map_rec_<ValueMap, ThisKey, ThisValue, Rest...> : public static_map_rec_<ValueMap, Rest...>
  {
    using This = static_map_rec_<ValueMap, ThisKey, ThisValue, Rest...>;
    using Base = static_map_rec_<ValueMap, Rest...>;
  
    friend struct lookup_item<This>;

    typename ValueMap::template type<ThisValue> value_;
  
  protected:
    template <typename Key>
    static constexpr auto find(Key k = Key{}) {
      if constexpr (std::is_same_v<ThisKey, Key>)
        return lookup_item<This>{};
      else
        return Base::find(k);
    }
    
  public:
    template <typename Action>
    void apply(const Action& f) {
      f(value_);
      Base::apply(f);
    }
    
    template <typename Return = void, typename Pred>
    Return* first(const Pred& p) {
      return p(value_)
        ? static_cast<Return*>(&value_)
        : Base::template first<Return, Pred>(p);
    }
  };
  

  template <typename... Args>
  class static_map_impl_ : public static_map_rec_<Args...>
  {
  public:
    // template <typename _Lookup, std::enable_if_t<std::is_base_of_v<lookup_item_tag, _Lookup>, int> = 0>
    // auto& operator[](_Lookup) {
    //   return _Lookup::get(this);
    // }

    // template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    // auto& operator[](Key k) {
    //   // using lookup = decltype(static_map<_Items...>::template find<_Key>());
    //   // static_assert(lookup{}, "Type key is not available");
    //   // return lookup::get(this);
    //   constexpr auto lookup = static_map_rec_<Args...>::find(k);
    //   static_assert(lookup, "Type key is not available");
    //   return lookup(this);
    // }

    template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    const auto& operator[](Key k) const {
      static constexpr auto lookup = static_map_rec_<Args...>::find(k);
      static_assert(lookup, "Type key is not available");
      return lookup(this);
    }

    template <typename Key, std::enable_if_t<!std::is_base_of_v<lookup_item_tag, Key>, int> = 0>
    auto& operator[](Key k) {
      static constexpr auto lookup = static_map_rec_<Args...>::find(k);
      static_assert(lookup, "Type key is not available");
      return lookup(this);
    }
  };


  template <typename... Rest>
  class static_map : public static_map_impl_<type_map::identity, Rest...> { };

  // template <typename... Rest>
  // class static_map<type_map::identity, Rest...> : public static_map_impl_<type_map::identity, Rest...> { };
  //
  // template <typename... Rest>
  // class static_map<type_map::unique_ptr, Rest...> : public static_map_impl_<type_map::unique_ptr, Rest...> { };
}