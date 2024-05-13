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

#include "ScalarConvert.hpp"

namespace dpl::qt::property_editor
{
  template <typename Scalar>
  struct DefaultFormat
  {
    static std::string Format(const Scalar& x) {
      return ScalarConvert<Scalar>::ToQString(x).toStdString();
    }
  };

  template <typename Scalar, int N>
  struct DefaultFormat<vector_n<Scalar, N>>
  {
    template <int i>
    static std::string Format(const Scalar& x) {
      static_assert(0 <= i && i < N, "Index is outside of the bounds");
      return DefaultFormat<Scalar>::Format(x);
    }
  };

  template <typename Derived, typename ValueType>
  struct BaseItem : DefaultFormat<ValueType>
  {
    using Type = ValueType;

    constexpr bool IsVisible() { return true; }
    constexpr bool Calculated() { return false; }
    constexpr bool IsReadOnly() { return false; }
    constexpr bool IsActive() { return false; }
    constexpr auto Tooltip() { return ""; }

    void SetActive() {}

    void Reset() {}

    template <typename T>
    void Set_(T&&) {}

    template <typename T>
    void Set(T&& v) {
      static_cast<Derived*>(this)->Set_(std::forward<T>(v));
    }

    Type Get() { return static_cast<Derived*>(this)->Get_(); }
  };


  template <typename ValueType>
  struct ItemFunctor : DefaultFormat<ValueType>
  {
    using Type = ValueType;

  private:
    std::string name;
    std::function<Type()> get;
    std::function<void(const Type&)> set;

  public:
    constexpr bool IsVisible() { return true; }
    constexpr bool Calculated() { return false; }
    bool IsReadOnly() { return !set; }
    constexpr bool IsActive() { return false; }
    constexpr auto Tooltip() { return ""; }

    void SetActive() {}
    void Reset() {}

    auto Name() {
      return name;
    }

    Type Get() {
      return get();
    }

    void Set(const Type& x) {
      set(x);
    }

    ItemFunctor(std::string name, const auto& get)
      : name(std::move(name)), get(get) {}

    ItemFunctor(std::string name, const auto& get, const auto& set)
      : name(std::move(name)), get(get), set(set) {}
  };
}
