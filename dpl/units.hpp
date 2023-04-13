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

#include <ratio>

namespace dpl::units
{
  template <typename Quantity>
  struct quantity_traits;
  

  inline static constexpr auto milli = std::milli{};
  inline static constexpr auto centi = std::centi{};
  

  constexpr auto pow(double value, int exponent) {
    auto result = 1.0;
    for (auto i = 0; i < exponent; ++i)
      result *= value;
    return result;
  }


  struct _no_coef_t {};
  
  template <typename Unit>
  inline constexpr auto coef = _no_coef_t{};

  template <typename Unit>
  inline constexpr bool has_coef = !std::is_same_v<decltype(coef<Unit>), _no_coef_t>;
  
  

  struct _unit_tag {};
  
  template<typename Unit, typename Quantity>
  struct _unit_base : _unit_tag
  {
  private:
    template<typename U>
    constexpr double parse(U u) {
      if constexpr (std::is_arithmetic_v<U>)
        return u;
      else {
        static_assert(std::is_same_v<quantity, typename U::quantity>, "Different quantities.");

        using SI = typename unit::quantity::SI;
        
        if constexpr (std::is_same_v<unit, SI> && has_coef<U>)
          return u.value*coef<U>;
        else if constexpr (std::is_same_v<U, SI> && has_coef<unit>)
          return u.value/coef<unit>;
        else if constexpr (has_coef<U> && has_coef<unit>)
          return u.value*coef<U>/coef<unit>;
        else
          static_assert(always_false<U>, "Cannot be converted.");
      }
    }
    
  public:
    using unit = Unit;
    using quantity = Quantity;

    double value = 1.0;

    constexpr operator double() const { return value; }

    const double& operator*() const { return value; }
    double& operator*() { return value; }
    
    template<typename U = double>
    constexpr _unit_base(U u = 1.0) {
      value = parse(u);
    }
    
    template<typename U, std::enable_if_t<std::is_arithmetic_v<U> || std::is_base_of_v<U, _unit_tag>, int> = 0>
    auto& operator=(U u) {
      value = parse(u);
      return static_cast<unit&>(*this); // NOLINT (misc-unconventional-assign-operator)
    }
    
    // auto& operator=(double v) {
    //   value = v;
    //   return static_cast<unit&>(*this); // NOLINT (misc-unconventional-assign-operator)
    // }

    // template <typename T>
    // auto& operator=(const T& t) {
    //   value = static_cast<double>(t);
    //   return static_cast<unit&>(*this); // NOLINT (misc-unconventional-assign-operator)
    // }
    
    constexpr unit operator*(double x) const {
      return {value*x};
    }

    constexpr unit operator/(double x) const {
      return {value/x};
    }
  };

  template<typename Unit, typename Quantity>
  constexpr auto operator*(double l, const _unit_base<Unit, Quantity>& r) {
    return Unit{l*r.value};
  }
  
  template<typename Unit, typename Quantity>
  constexpr auto operator/(double l, const _unit_base<Unit, Quantity>& r) {
    return Unit{l/r.value};
  }

#ifdef INCLUDE_NLOHMANN_JSON_HPP_
  template <typename Unit, typename Quantity>
  void to_json(nlohmann::json& j, const _unit_base<Unit, Quantity>& u) {
    j = u.value;
  }

  template <typename Unit, typename Quantity>
  void from_json(const nlohmann::json& j, _unit_base<Unit, Quantity>& u) {
    u.value = j;
  }
#endif



  

  // NOLINTNEXTLINE (cppcoreguidelines-macro-usage)
  #define DPL_UNIT(unit, quantity) \
    struct unit : _unit_base<unit, quantity> { using _unit_base<unit, quantity>::operator=; }; 

  template<typename UnitSI, int Dim = 1>
  struct _quantity_base
  {
    using SI = UnitSI;
    static constexpr auto dim = Dim;
  };
  
  // -------
    
  struct length;

  DPL_UNIT(metre, length)
  DPL_UNIT(centimetre, length)
  DPL_UNIT(foot, length)
  DPL_UNIT(inch, length)
  
  struct length : _quantity_base<metre> {};

  // -------
  
  struct volume;

  DPL_UNIT(metre3, volume)
  DPL_UNIT(foot3, volume)
  DPL_UNIT(barrel, volume)

  struct volume : _quantity_base<metre3, 3> {};

  // -------
  
  struct pressure;

  DPL_UNIT(pascal, pressure)
  DPL_UNIT(megapascal, pressure)
  DPL_UNIT(bar, pressure)
  DPL_UNIT(psi, pressure)

  struct pressure : _quantity_base<pascal> {};

  // -------

  struct time;

  DPL_UNIT(second, time)
  DPL_UNIT(day, time)
  DPL_UNIT(year, time)
  
  struct time : _quantity_base<second> {};

  // -------

  struct viscosity;

  DPL_UNIT(pascal_second, viscosity)
  DPL_UNIT(poise, viscosity)
  DPL_UNIT(centipoise, viscosity)

  struct viscosity : _quantity_base<pascal_second> {};

  // -------



  
  

  

  template <> inline constexpr double coef<foot> = 0.3048;
  template <> inline constexpr double coef<centimetre> = 0.01;
  template <> inline constexpr double coef<inch> = coef<foot>/12;
  
  template <> inline constexpr double coef<foot3> = pow(coef<foot>, 3);
  template <> inline constexpr double coef<barrel> = 0.158987;
  
  template <> inline constexpr double coef<day> = 86400;
  template <> inline constexpr double coef<year> = 365*coef<day>;

  template <> inline constexpr double coef<bar> = 1e5;
  template <> inline constexpr double coef<megapascal> = 1e6;
  template <> inline constexpr double coef<psi> = 6894.7572932;

  template <> inline constexpr double coef<poise> = 0.1;
  template <> inline constexpr double coef<centipoise> = 0.01*coef<poise>;



  

  template<std::intmax_t N, intmax_t D, typename U, std::enable_if_t<std::is_base_of_v<_unit_tag, U>, int> = 0>
  constexpr auto operator^(const std::ratio<N, D>& l, const U& r) {
    return r*pow(1.0*l.num/l.den, U::quantity::dim);
  }

  template<typename Stream, typename U, std::enable_if_t<std::is_base_of_v<_unit_tag, U>, int> = 0>
  auto& operator<<(Stream& s, const U& r) {
    s << r.value;
    return s;
  }

  
  
       
  
  
  

  
}
