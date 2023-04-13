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
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <cmath>


namespace dpl
{
  template <typename... T>
  constexpr bool always_false = false;
  
  template <typename... Args_>
  struct has_any : std::bool_constant<sizeof...(Args_) != 0> {};

  template <typename... Args_>
  inline constexpr bool has_any_v = has_any<Args_...>::value;

  template <int i_>
  using int_const = std::integral_constant<int, i_>;


  
  template <typename V_, typename ...Rest_>
  struct are_arithmetic : std::bool_constant<std::is_arithmetic_v<V_> && are_arithmetic<Rest_...>::value> {};

  template <typename V_>
  struct are_arithmetic<V_> : std::is_arithmetic<V_> {};

  
  template <typename U_, typename V_, typename ...Rest_>
  struct are_assignable : std::bool_constant<std::is_assignable_v<U_, V_> && are_assignable<U_, Rest_...>::value> {};

  template <typename U_, typename V_>
  struct are_assignable<U_, V_> : std::is_assignable<U_, V_> {};


  

  

  // template <typename... Args_>
  // inline constexpr bool are_arithmetic_v = are_arithmetic<Args_...>::value;


  // std::bool_constant<std::is_arithmetic_v<V_> && are_arithmetic<Rest_...>::value>
  // std::bool_constant<std::is_arithmetic_v<V_> && are_arithmetic<Rest_...>::value> {};

  


  template <int First, int Last, typename Lambda>
  constexpr void sfor(const Lambda& f) {
    if constexpr (First < Last) {
      if constexpr (std::is_invocable_v<Lambda, int_const<First>>)
        f(int_const<First>{});
      else
        f();

      dpl::sfor<First + 1, Last>(f);
    }
  }

  template <int Count, typename Lambda>
  constexpr void sfor(const Lambda& f) {
    dpl::sfor<0, Count>(f);
  }

  // template<int Count>
  // struct static_dim
  // {
  //   template<int Idx>
  //   static constexpr auto next() {
  //     if constexpr (Idx < Count - 1)
  //       return Idx + 1;
  //     else
  //       return 0;
  //   }
  //
  //   template<int Idx>
  //   static constexpr auto next(std::integral_constant<int, Idx>) {
  //     return std::integral_constant<int, next<Idx>()>{};
  //   }  
  // };


  
  struct extrapolant
  {
    struct linear {};
    struct flat {};
  };

  struct lerp_base
  {
    struct linear {};
    struct log10 {};
  };
  
  
  
  template<typename Base = lerp_base::linear, typename Extrapolant = extrapolant::linear, typename T>
  constexpr auto lerp(const T& a, const T& b, double t) {
    if constexpr (std::is_same_v<Extrapolant, extrapolant::linear> && std::is_same_v<Base, lerp_base::linear>)
      return a + (b - a)*t;
    else if constexpr (std::is_same_v<Base, lerp_base::log10>)
      return std::pow(10., dpl::lerp<lerp_base::linear, Extrapolant>(std::log10(a), std::log10(b), t));
    else if constexpr (std::is_same_v<Extrapolant, extrapolant::flat>)
      return
        t < 0 ? a :
        t > 1 ? b :
        dpl::lerp<Base, extrapolant::linear>(a, b, t);
    else
      static_assert(always_false<T>, "wrong extrapolant or base");
  }

  template<typename Base = lerp_base::linear, typename Extrapolant = extrapolant::linear, typename T>
  constexpr auto lerp(T* ptr, double t) {
    return dpl::lerp<Base, Extrapolant>(ptr[0], ptr[1], t);
  }

  template<typename Extrapolant = extrapolant::linear, typename T>
  constexpr auto lerp(bool log10, const T& a, const T& b, double t) {
    return log10
      ? dpl::lerp<lerp_base::log10, Extrapolant>(a, b, t)
      : dpl::lerp<lerp_base::linear, Extrapolant>(a, b, t);
  }

  template<typename Extrapolant = extrapolant::linear, typename T>
  constexpr auto lerp(bool log10, T* ptr, double t) {
    return dpl::lerp<Extrapolant>(log10, ptr[0], ptr[1], t);
  }
  


  template<char thousands_sep_ = ',', typename T_>
  std::string ft_format(const T_& value, int precision = 0) {
    struct separate_thousands : std::numpunct<char> {
      char_type do_thousands_sep() const override { return thousands_sep_; }  // separate with commas
      string_type do_grouping() const override { return "\3"; } // groups of 3 digit
    };

    std::ostringstream ss;
    ss.imbue(std::locale(ss.getloc(), new separate_thousands));
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
  }
  



  // template <typename Func>
  // void apply(const Func& f) {
  //   f(value_);
  //   Base::apply(f);
  // }

  // template <typename Action>
  // void apply(const Action&) {}

  // template<typename Action, typename Item, typename... Rest>
  // void apply(const Action& a, Item&& item, Rest&&... args) {
  //   a(std::forward<Item>(item));
  //   dpl::apply(a, std::forward<Rest>(args)...);
  // }
  
  // inline void AddItem(QBoxLayout* ) { }
  //
  // template<typename... Rest>
  // void AddItem(QBoxLayout* l, QWidget& widget, Rest&&... args) {
  //   l->addWidget(&widget);
  //   QLayoutBuilder::AddItem(l, std::forward<Rest>(args)...);
  // }
}





