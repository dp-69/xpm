/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2024
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

#include "units.hpp"
#include "static_vector.hpp"

#include <fmt/core.h>

#include <filesystem>

template <>
struct fmt::formatter<std::filesystem::path> : formatter<std::string_view>
{
  template <typename FormatContext>
  auto format(const std::filesystem::path& path, FormatContext& ctx) {
    return formatter<std::string_view>::format(path.string(), ctx);
  }
};

template <typename T, typename Tag>
struct fmt::formatter<dpl::strong_integer<T, Tag>> : formatter<T>
{
  template <typename FormatContext>
  auto format(dpl::strong_integer<T, Tag> si, FormatContext& ctx) {
    return formatter<T>::format(*si, ctx);
  }
};

template <typename T>
struct fmt::formatter<dpl::vector_n<T, 3>>
{
  template<typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(dpl::vector_n<T, 3> si, FormatContext& ctx) {
    return fmt::format_to(ctx.out(), "({}, {}, {})", si.x(), si.y(), si.z());
  }
};

template <typename Quantity, typename Tag>
struct fmt::formatter<dpl::units::def_unit<Quantity, Tag>> // : formatter<double>
{
  std::string fmt;

  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
      std::string tmp;
      auto it = ctx.begin();
      while(it != ctx.end() && *it != '}')
        tmp += *it++;
      fmt = tmp;
      return it;
  }

  template <typename FormatContext>
  auto format(dpl::units::def_unit<dpl::units::digital_size, Tag> si, FormatContext& ctx) {
    if (fmt.empty())
      return fmt::format_to(ctx.out(), "{:.0f}", *si);

    return fmt::format_to(ctx.out(), fmt::runtime("{:"+fmt+"}"), *si);
  }
};