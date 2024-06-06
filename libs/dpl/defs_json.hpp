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
 * along with DPL. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "units.hpp"
// #include "static_vector.hpp"
#include "curve2d.hpp"

#include <nlohmann/json.hpp>

namespace dpl
{
  template <typename Tag>
  void to_json(nlohmann::json& j, const so_integer<Tag>& i) {
    j = i.value;
  }

  template <typename Tag>
  void from_json(const nlohmann::json& j, so_integer<Tag>& i) {
    i.value = j;
  }

  namespace units
  {
    template <typename Quantity, typename Tag>
    void to_json(nlohmann::json& j, const def_unit<Quantity, Tag>& unit) {
      j = unit.value;
    }

    template <typename Quantity, typename Tag>
    void from_json(const nlohmann::json& j, def_unit<Quantity, Tag>& unit) {
      unit.value = j;
    }
  }

  // ----------------------------------

  template <typename Derived, typename T, int n>
  void to_json(nlohmann::json& j, const _vector_oper<Derived, T, n>& v) {
    sfor<n>([&](auto i) { j.push_back(v[i]); });
  }

  template <typename Derived, typename T, int n>
  void from_json(const nlohmann::json& j, _vector_oper<Derived, T, n>& v) {
    if (j.is_array())
      sfor<n>([&](auto i) { v[i] = j[i]; });
    else { // if (j.is_number())
      sfor<n>([&](auto i) { v[i] = j; });
    }
  }

  template <typename T, int n>
  void parse(const nlohmann::json& j, std::vector<vector_n<T, n>>& v) {
    if (j.empty())
      v.clear();

    if (j[0].is_array()) {

    }
    else {
      auto size = j.size()/n;
      v.resize(size);
      for (size_t i = 0; i < size; ++i)
        dpl::sfor<n>([&](auto d) {
          v[i][d] = j[i*n + d];
        });
        // for (auto d = 0; d < n; ++d)
        //   v[i][d] = j[i*n + d];
    }
  }

  inline void parse(const nlohmann::json& j, curve2d& v) {
    parse(j, v.storage);
  }

  // ---------------------------------


  inline nlohmann::json parse_path(const nlohmann::json& j, bool ignore_comments = false) {
    return j.is_string() ? nlohmann::json::parse(std::ifstream{j.get<std::string>()}, nullptr, true, ignore_comments) : j;
  }

  template<typename T>
  concept json_parsable = requires (const nlohmann::json& j, T& v) { parse(j, v); };

  template<typename T>
  concept not_json_parsable = !json_parsable<T>;

  // NOLINTBEGIN(bugprone-empty-catch)

  template <json_parsable T>
  void try_parse(const nlohmann::json& j, T& v) {
    try {
      parse(j, v);
    } catch (...) {}
  }

  template <json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, T& v) {
    try {
      parse(j.at(k), v);
    } catch (...) {}
  }

  template <json_parsable T, typename Key, int n>
  void try_parse(const nlohmann::json& j, const Key& k, std::array<T, n>& v) {
    try {
      const auto& jk = j.at(k);

      for (int i = 0; i < n; ++i)
        parse(jk[i], v[i]);
    } catch (...) {}
  }

  // template <not_json_parsable T, typename Key, int n>
  // void try_parse(const nlohmann::json& j, const Key& k, std::array<std::optional<T>, n>& v) {
  //   try {
  //     const auto& jk = j.at(k);
  //
  //     for (int i = 0; i < n; ++i)
  //       if (const auto& x = jk[i]; x.is_null())
  //         v[i].reset();
  //       else
  //         v[i] = x;
  //   } catch (...) {}
  // }

  template <json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, T& v, bool ignore_comments) {
    try {
      parse(parse_path(j.at(k), ignore_comments), v);
    } catch (...) {}
  }

  template <not_json_parsable T>
  void try_parse(const nlohmann::json& j, T& v) {
    try { 
      v = j;
    } catch (...) {}
  }

  template <not_json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, T& v) {
    try {
      v = j.at(k);
    } catch (...) {}
  }

  template <not_json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, std::optional<T>& v) {
    try {
      if (const nlohmann::json& jj = j.at(k); jj.is_null())
        v.reset();
      else
        v = jj;
    } catch (...) {}
  }

  template <not_json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, std::optional<T>& v, bool ignore_comments) {
    try {
      if (const nlohmann::json& jj = parse_path(j.at(k), ignore_comments); jj.is_null())
        v.reset();
      else
        v = jj;
    } catch (...) {}
  }

  template <not_json_parsable T, typename Key, int n>
  void try_parse(const nlohmann::json& j, const Key& k, std::array<std::optional<T>, n>& v) {
    try {
      const auto& jk = j.at(k);

      for (int i = 0; i < n; ++i)
        if (const auto& x = jk[i]; x.is_null())
          v[i].reset();
        else
          v[i] = x;
    } catch (...) {}
  }

  template <typename GetType, not_json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, T& v) {
    try {
      v = j.at(k).template get<GetType>();
    } catch (...) {}
  }

  template <typename GetType, not_json_parsable T, typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, std::optional<T>& v) {
    try {
      if (const nlohmann::json& jj = j.at(k); jj.is_null())
        v.reset();
      else
        v = jj.get<GetType>();
    } catch (...) {}
  }

  template <typename Key>
  void try_parse(const nlohmann::json& j, const Key& k, nlohmann::json& v, bool ignore_comments) {
    try {
      v = parse_path(j.at(k), ignore_comments);
    } catch (...) {}
  }

  // NOLINTEND(bugprone-empty-catch)
}