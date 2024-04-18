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

#include "general.hpp"

#include <concepts>

namespace dpl::units
{
  template <typename>
  inline constexpr auto coef = nullptr;

  template <typename T>
  inline constexpr bool has_coef = std::is_arithmetic_v<decltype(coef<T>)>;

  template<typename Quantity, typename Tag>
  struct def_unit {
    template <typename U> requires std::is_arithmetic_v<U>
    static constexpr auto parse(const U u) {
      return u;
    }

    template <typename U>                                                            
    static constexpr auto parse(const def_unit<Quantity, U> u) {
      using standard_tag = typename Quantity::standard;                                      
                                                                                     
      if constexpr (std::is_same_v<Tag, standard_tag> && has_coef<U>)                         
        return u.value*coef<U>;                                                     
      else if constexpr (std::is_same_v<U, standard_tag> && has_coef<Tag>)                    
        return u.value/coef<Tag>;                                                  
      else if constexpr (has_coef<U> && has_coef<Tag>)                              
        return u.value*coef<U>/coef<Tag>;                                          
      else                                                                           
        static_assert(always_false<U>, "Cannot be converted.");                      
    }

    double value;

    auto& operator*() const { return value; }
    auto& operator*() { return value; }

    constexpr def_unit() : value{1.0} {}                                                                            

    // constexpr def_unit(const auto u)  {                                              
    //   value = parse(u);                                                    
    // }    

    template <typename U> requires std::is_arithmetic_v<U>
    constexpr def_unit(const U u) {                                              
      value = parse(u);                                                    
    }
    
    template <typename U>
    constexpr def_unit(const def_unit<Quantity, U> u) {                                              
      value = parse(u);                                                    
    }

    template <typename U> requires std::is_arithmetic_v<U>                       
    auto& operator=(const U u) {                                                 
      value = parse(u);                                                    
      return *this;                                                              
    }                                                                            
                                                                                 
    template <typename U>                                                        
    auto& operator=(const def_unit<Quantity, U> u) {                                 
      value = parse(u);                                                    
      return *this;                                                              
    }
  };

  

  template<typename Quantity, typename Tag>
  inline constexpr auto coef<def_unit<Quantity, Tag>> = coef<Tag>;

  // template <typename>
  // struct quantity_dim { static constexpr int value = 1; };

  // #define DEF_DIM(quantity, dim) \
  //   template <> struct quantity_dim<quantity> { static constexpr int value = dim; };

  // NOLINTBEGIN(bugprone-macro-parentheses)
  #define DEF_SI(quantity, unit)                                              \
    namespace helper { struct unit##_tag {}; }                                \
    struct quantity { using standard = helper::unit##_tag; };                 \
    using unit = def_unit<quantity, helper::unit##_tag>;


  #define DEF_NON_SI(quantity, unit, mult)                                    \
    namespace helper { struct unit##_tag {}; }                                \
    using unit = def_unit<quantity, helper::unit##_tag>;                      \
    template <> inline constexpr double coef<helper::unit##_tag> = mult; 
  // NOLINTEND(bugprone-macro-parentheses)

  DEF_SI(length, metre)
    DEF_NON_SI(length, centimetre, 0.01)
    DEF_NON_SI(length, foot, 0.3048)
    DEF_NON_SI(length, inch, coef<foot>/12)
  
  DEF_SI(volume, metre3)
    DEF_NON_SI(volume, foot3, pow(coef<foot>, 3))
    DEF_NON_SI(volume, barrel, 0.158987)

  DEF_SI(pressure, pascal)
    DEF_NON_SI(pressure, megapascal, 1e6)
    DEF_NON_SI(pressure, bar, 1e5)
    DEF_NON_SI(pressure, psi, 6894.7572932)

  DEF_SI(time, second)
    DEF_NON_SI(time, day, 86400)
    DEF_NON_SI(time, year, 365*coef<day>)
  
  DEF_SI(viscosity, pascal_second)
    DEF_NON_SI(viscosity, poise, 0.1)
    DEF_NON_SI(viscosity, centipoise, 0.01*coef<poise>)

  DEF_SI(area, metre2)
    DEF_NON_SI(area, darcy, 9.869233e-13)
    DEF_NON_SI(area, millidarcy, 1e-3*coef<darcy>)

  using permeability = area;

  DEF_SI(digital_size, byte)
    DEF_NON_SI(digital_size, megabyte, pow(2, 20))
    DEF_NON_SI(digital_size, gigabyte, pow(2, 30))

  // -------

  struct pore_volume
  {
    double value;

    pore_volume(const double value) : value(value) {}

    auto& operator*() const { return value; }
    auto& operator*() { return value; }
  };


  template <typename Stream, typename Quantity, typename Tag>
  auto& operator<<(Stream& stream, const def_unit<Quantity, Tag>& unit) {
    stream << unit.value;
    return stream;
  }

  
  #undef DEF_SI
  #undef DEF_UNIT
}