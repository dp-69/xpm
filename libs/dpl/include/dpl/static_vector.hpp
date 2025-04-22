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

#include "general.hpp"

#include <ostream>
#include <span>
#include <algorithm>

namespace dpl
{
  namespace detail
  {
    struct op
    {
      static constexpr auto add  = [](const auto& l, const auto& r) { return l + r; };
      static constexpr auto diff = [](const auto& l, const auto& r) { return l - r; };
      static constexpr auto mult = [](const auto& l, const auto& r) { return l*r; };
      static constexpr auto div  = [](const auto& l, const auto& r) { return l/r; };
      
      struct cross
      {
        static constexpr std::integral_constant<int, 0> _0{};
        static constexpr std::integral_constant<int, 1> _1{};
        static constexpr std::integral_constant<int, 2> _2{};
        
        static constexpr auto x = [](const auto& a, const auto& b) { return a[_1]*b[_2] - a[_2]*b[_1]; };
        static constexpr auto y = [](const auto& a, const auto& b) { return a[_2]*b[_0] - a[_0]*b[_2]; };
        static constexpr auto z = [](const auto& a, const auto& b) { return a[_0]*b[_1] - a[_1]*b[_0]; };
      };
    };
  }
  
  
  template <typename T_, int n_>
  class vector_n;
  
  template <typename Derived_, typename T_, int n_>
  class _vector_oper
  {
    static constexpr std::integral_constant<int, 0> _0{};
    
    template <typename U_, typename Func_>
    constexpr auto _binary_op_impl_ptr(U_* r, const Func_& f) const {
      vector_n<decltype(f(T_{}, U_{})), n_> v;
      sfor<n_>([&](auto i) { v[i] = f((*this)[i], r[i]); });
      return v;
    }

    template <typename U_, typename Func_>
    constexpr auto _binary_op_impl(const U_& r, const Func_& f) const {
      if constexpr (std::is_array_v<U_> || std::is_pointer_v<U_>)
        return _binary_op_impl_ptr(r, f);
      else if constexpr (std::is_arithmetic_v<U_>)
        return vector_n<decltype(f(T_{}, U_{})), n_>{f, *this, r};
      else
        return vector_n<decltype(f(T_{}, r[_0])), n_>{f, *this, r};
    }

  public:
    using ScalarType = T_;
    
    // template <int i>           auto& get()       { return static_cast<      Derived_*>(this)->template _get<i>(); }
    // template <int i> constexpr auto& get() const { return static_cast<const Derived_*>(this)->template _get<i>(); }
    
    template <int i> auto& operator[](std::integral_constant<int, i>) {
      return static_cast<Derived_*>(this)->template get<i>();
    }
    
    template <int i> constexpr auto& operator[](std::integral_constant<int, i>) const {
      return static_cast<const Derived_*>(this)->template get<i>();
    }

    template <typename U = void>
    constexpr auto prod() const /*-> std::conditional_t<std::is_same_v<U, void>, decltype(operation::mult(T_{}, T_{})), U>*/ {
      return dpl::tail_aggregate<n_, U>(detail::op::mult, *this);
    }

    constexpr auto operator+(const auto& r) const { return _binary_op_impl(r, detail::op::add); }
    constexpr auto operator-(const auto& r) const { return _binary_op_impl(r, detail::op::diff); }
    constexpr auto operator*(const auto& r) const { return _binary_op_impl(r, detail::op::mult); }
    constexpr auto operator/(const auto& r) const { return _binary_op_impl(r, detail::op::div); }
  };

  template <typename T, int n>
  constexpr auto cross(const vector_n<T, n>& lhs, const vector_n<T, n>& rhs) {
    return vector_n<T, n>{detail::op::cross{}, lhs, rhs};
  }

  template <typename T>
  constexpr int non_zero_idx(const vector_n<T, 3>& v) {
    if (v.x) return 0;
    if (v.y) return 1;
    if (v.z) return 2;
    return -1;
  }

  template <typename Lhs, typename Derived, typename T, int n> requires (std::is_arithmetic_v<Lhs>)
  constexpr auto operator*(Lhs lhs, const _vector_oper<Derived, T, n>& rhs) {
    return rhs*lhs;
  }

  template <typename Derived, int n>
  constexpr auto operator/(double l, const _vector_oper<Derived, double, n>& r) { // TODO
    vector_n<double, n> v;
    dpl::sfor<n>([&](auto i) { v[i] = l/r[i]; });
    return v;
  }
  
  template <typename Derived, typename T, int n>
  std::ostream& operator<<(std::ostream& os, const _vector_oper<Derived, T, n>& v) {
              os << v[std::integral_constant<int, 0>{}];
    if constexpr (n > 1)
      os << ", " << v[std::integral_constant<int, 1>{}];
    if constexpr (n > 2)
      os << ", " << v[std::integral_constant<int, 2>{}];

    return os;
  }
  

  

  
  template <typename T>
  class vector_n<T, 2> : public _vector_oper<vector_n<T, 2>, T, 2>
  {
  public:
    T x, y;

    template <int i>
    constexpr auto& get() {
           if constexpr (i == 0) return x;
      else if constexpr (i == 1) return y;
    }

    template <int i>
    constexpr auto& get() const {
           if constexpr (i == 0) return x;
      else if constexpr (i == 1) return y;
    }

  private:
    static constexpr std::integral_constant<int, 0> _0{};
    static constexpr std::integral_constant<int, 1> _1{};
    static constexpr std::integral_constant<int, 2> _2{};
    
  public:
    operator T*() { return &x; }
    operator const T*() const { return &x; }

    //       T* ptr()       { return &x; }
    // const T* ptr() const { return &x; }
    
    vector_n() = default;

    template <typename U> requires (std::is_same_v<U, T> || std::is_arithmetic_v<U>)
    constexpr vector_n(const U val) {
      x = val;
      y = val;
    }

    template <typename X, typename Y> requires (std::is_assignable_v<T&, X> && std::is_assignable_v<T&, Y>)
    constexpr vector_n(const X x0, const Y y0) {
      x = x0;
      y = y0;
    }
    
    template <typename Func> requires std::is_invocable_v<Func, decltype(_0)>
    explicit constexpr vector_n(const Func& f) : x{f(_0)}, y{f(_1)} {}

    template <typename U, typename V, typename Func>
      requires std::is_invocable_v<Func, typename U::ScalarType, typename V::ScalarType>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, 2>& l,
      const _vector_oper<V, typename V::ScalarType, 2>& r
    ) : x{f(l[_0], r[_0])}, y{f(l[_1], r[_1])} {}

    template <typename U, typename V, typename Func> requires std::is_arithmetic_v<V>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, 2>& lhs,
      const V rhs) : x{f(lhs[_0], rhs)}, y{f(lhs[_1], rhs)} {}

    template <typename Derived, typename U>
    constexpr vector_n(const _vector_oper<Derived, U, 2>& v)
      : x{static_cast<T>(v[_0])}, y{static_cast<T>(v[_1])} {}
    
    template <int i, typename U = int>
    explicit constexpr vector_n(std::integral_constant<int, i>, const U v = 1) {
           if constexpr (i == 0) { x = v; y = 0; }
      else if constexpr (i == 1) { x = 0; y = v; }
    }
    
    vector_n(auto* ptr) : x{ptr[0]}, y{ptr[1]} {}
    vector_n(const auto* ptr) : x{ptr[0]}, y{ptr[1]} {}

    auto length() const { return std::sqrt(x*x + y*y); }

    auto normalise() const {
      auto l = length();
      return vector_n{x/l, y/l};
    }

    constexpr auto slope() const { return y/x; }

    constexpr auto min() const { return std::min(x, y); }
    constexpr auto sum() const { return x + y; }
  };

  template <typename T>
  class vector_n<T, 3> : public _vector_oper<vector_n<T, 3>, T, 3>
  {
  public:
    T x, y, z;

    template <int i>
    constexpr auto& get() {
           if constexpr (i == 0) return x;
      else if constexpr (i == 1) return y;
      else if constexpr (i == 2) return z;
    }

    template <int i>
    constexpr auto& get() const {
           if constexpr (i == 0) return x;
      else if constexpr (i == 1) return y;
      else if constexpr (i == 2) return z;
    }
    
  private:
    static constexpr std::integral_constant<int, 0> _0{};
    static constexpr std::integral_constant<int, 1> _1{};
    static constexpr std::integral_constant<int, 2> _2{};


  public:
    operator T*() { return &x; }
    operator const T*() const { return &x; }
    
    //       T* ptr()       { return &x; }
    // const T* ptr() const { return &x; }
    
    vector_n() = default;

    template <typename U> requires (std::is_same_v<U, T> || std::is_arithmetic_v<U>)
    constexpr vector_n(const U val) {
      x = val;
      y = val;
      z = val;
    }

    template <typename X, typename Y, typename Z> requires (
      std::is_assignable_v<T&, X> &&
      std::is_assignable_v<T&, Y> &&
      std::is_assignable_v<T&, Z>)
    constexpr vector_n(const X x0, const Y y0, const Z z0) {
      x = x0;
      y = y0;
      z = z0;
    }

    // template <typename ...Args> requires (sizeof...(Args) == 3 && are_assignable_v<T&, Args...>) // TODO
    // constexpr vector_n(Args... args) {
    //   _assign<Args...>(args...);
    // }

    template <typename Func> requires std::is_invocable_v<Func, decltype(_0)>
    explicit constexpr vector_n(const Func& f) : x{f(_0)}, y{f(_1)}, z{f(_2)} {}

    template <typename U, typename V, typename Func>
      requires std::is_invocable_v<Func, typename U::ScalarType, typename V::ScalarType>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, 3>& l,
      const _vector_oper<V, typename V::ScalarType, 3>& r
    ) : x{f(l[_0], r[_0])}, y{f(l[_1], r[_1])}, z{f(l[_2], r[_2])} {}

    template <typename U, typename V, typename Func> requires std::is_arithmetic_v<V>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, 3>& lhs,
      const V rhs) : x{f(lhs[_0], rhs)}, y{f(lhs[_1], rhs)}, z{f(lhs[_2], rhs)} {}

    template <typename Derived, typename U>
    constexpr vector_n(const _vector_oper<Derived, U, 3>& v)
      : x{static_cast<T>(v[_0])}, y{static_cast<T>(v[_1])}, z{static_cast<T>(v[_2])} {}

    template <typename U, typename V>
    constexpr vector_n(detail::op::cross, const U& l, const V& r)
      : x{detail::op::cross::x(l, r)},
        y{detail::op::cross::y(l, r)},
        z{detail::op::cross::z(l, r)} {}

    template <int i, typename U = int>
    explicit constexpr vector_n(std::integral_constant<int, i>, const U v = 1) {
           if constexpr (i == 0) { x = v; y = 0; z = 0; }
      else if constexpr (i == 1) { x = 0; y = v; z = 0; }
      else if constexpr (i == 2) { x = 0; y = 0; z = v; }
    }
    
    vector_n(auto* ptr) : x{ptr[0]}, y{ptr[1]}, z{ptr[1]} {}
    vector_n(const auto* ptr) : x{ptr[0]}, y{ptr[1]}, z{ptr[1]} {}

    template <typename Rhs = vector_n>
    constexpr auto dot(const Rhs& rhs) const { return x*rhs.x + y*rhs.y + z*rhs.z; }

    T length() const {
      if constexpr (std::is_arithmetic_v<T>)
        return std::sqrt(x*x + y*y + z*z);
      else
        return sqrt(x*x + y*y + z*z);
    }

    auto normalise() const {
      auto l = length();
      return vector_n{x/l, y/l, z/l};
    }

    constexpr auto min() const { return std::min(std::min(x, y), z); }
    constexpr auto max() const { return std::max(std::max(x, y), z); }
    constexpr auto sum() const { return x + y + z; }
  };
  
  
  template <typename T, int n>
  class vector_n_map : public _vector_oper<vector_n_map<T, n>, T, n>
  {
    T* ptr_;

    

    // friend class _vector_oper<vector_n_map, T, n>;

  public:
    template <int i>       auto& get()       { return ptr_[i]; }
    template <int i> const auto& get() const { return ptr_[i]; }
    
    operator T*() { return ptr_; }
    operator const T*() const { return ptr_; }

    vector_n_map(T* ptr) : ptr_(ptr) {}
  };
   

  

  
  using vector2i = vector_n<int, 2>;
  using vector2f = vector_n<float, 2>;
  using vector2d = vector_n<double, 2>;
  using vector2i_map = vector_n_map<int, 2>;
  using vector3i = vector_n<int, 3>;
  using vector3f = vector_n<float, 3>;
  using vector3d = vector_n<double, 3>;
  using vector3i_map = vector_n_map<int, 3>;
  using vector3d_map = vector_n_map<double, 3>;

  // template <size_t i, typename T, int n>
  //   decltype(auto) get(const dpl::vector_n<T, n>& v) {
  //          if constexpr (i == 0) return v.x;
  //     else if constexpr (i == 1) return v.y;
  //     else if constexpr (i == 2) return v.z;
  //   }
  //
  //   template <size_t i, typename T, int n>
  //   decltype(auto) get(dpl::vector_n<T, n>& v) {
  //          if constexpr (i == 0) return v.x;
  //     else if constexpr (i == 1) return v.y;
  //     else if constexpr (i == 2) return v.z;
  //   }
  //
  //   template <size_t i, typename T, int n>
  //   decltype(auto) get(dpl::vector_n<T, n>&& v) {
  //          if constexpr (i == 0) return v.x;
  //     else if constexpr (i == 1) return v.y;
  //     else if constexpr (i == 2) return v.z;
  //   }

  template <typename T>
  auto solve(const vector_n<T, 2>& p0, const vector_n<T, 2>& p1, const auto arg) {
    return p0.y + (p1 - p0).slope()*(arg - p0.x);
  }

  template <typename T>
  auto solve(const vector_n<T, 2>& p0, const vector_n<T, 2>& p1, const auto arg, loga_t, linear_t) {
    /* linear-log scale
     *
     * auto log_p0_y = std::log10(p0.y);
     * auto slope = (std::log10(p1.y) - log_p0_y)/(p1.x - p0.x);
     * return std::pow(10, log_p0_y + slope*(arg - p0.x));
     *
     *
     */

    auto slope = (p1.y - p0.y)/(std::log10(p1.x) - std::log10(p0.x));
    return p0.y + slope*(std::log10(arg) - std::log10(p0.x));

    // return p0.y + (p1 - p0).slope()*(arg - p0.x);
  }

  template <typename T, size_t Extent = std::dynamic_extent>
  auto solve(std::span<const vector_n<T, 2>, Extent> curve, const auto arg, linear_t = linear) {
    if (arg <= curve.front().x)
      return solve(curve[0], curve[1], arg);

    if (arg > curve.back().x)
      return solve(curve[curve.size() - 2], curve[curve.size() - 1], arg);

    auto iter = std::ranges::lower_bound(curve, arg, {}, [](const vector_n<T, 2>& p) { return p.x; });
    return solve(*(iter - 1), *iter, arg);
  }

  template <typename T, size_t Extent = std::dynamic_extent>
  auto solve(std::span<const vector_n<T, 2>, Extent> curve, const auto arg, flat_t) {
    if (arg <= curve.front().x)
      return curve.front().y;

    if (arg > curve.back().x)
      return curve.back().y;

    auto iter = std::ranges::lower_bound(curve, arg, {}, [](const vector_n<T, 2>& p) { return p.x; });
    return solve(*(iter - 1), *iter, arg);
  }

  template <typename T, size_t Extent = std::dynamic_extent>
  auto solve(std::span<const vector_n<T, 2>, Extent> curve, const auto arg, flat_t, loga_t, linear_t) {
    if (arg <= curve.front().x)
      return curve.front().y;

    if (arg > curve.back().x)
      return curve.back().y;

    auto iter = std::ranges::lower_bound(curve, arg, {}, [](const vector_n<T, 2>& p) { return p.x; });
    return solve(*(iter - 1), *iter, arg, loga, linear);
  }

  // template <typename T>
  // auto solve(const std::vector<vector_n<T, 2>>& curve, const auto arg, auto extrapolant) {
  //   return solve(std::span{curve}, arg, extrapolant);
  // }
}


// NOLINTBEGIN(cert-dcl58-cpp)
namespace std
{
  template <typename Type, int n>
  struct tuple_size<dpl::vector_n<Type, n>> : integral_constant<size_t, n> {};
  
  template <size_t i, typename T, int n>
  struct tuple_element<i, dpl::vector_n<T, n>> { using type = T; };
  
  template <typename T, int n>
  struct tuple_size<dpl::vector_n_map<T, n>> : integral_constant<size_t, n> {};    
  
  template <size_t i, typename T, int n>
  struct tuple_element<i, dpl::vector_n_map<T, n>> { using type = T; };

  /* 
   * fix to support std::ranges::views::elements_view in c++20
   *
   * 
   * template< class T, std::size_t N >
   * concept  =
   *     requires(T t) {
   *         typename std::tuple_size<T>::type;
   *         requires N < std::tuple_size_v<T>;
   *         typename std::tuple_element_t<N, T>;
   *         { std::get<N>(t) } -> std::convertible_to<const std::tuple_element_t<N, T>&>;
   *     };
   *
   */
  template <size_t i, typename T, int n>
  const T& get(const dpl::vector_n<T, n>& v) {
         if constexpr (i == 0) return v.x;
    else if constexpr (i == 1) return v.y;
    else if constexpr (i == 2) return v.z;
  }

  // template <size_t i, typename T, int n>
  // T& get(dpl::vector_n<T, n>& v) {
  //        if constexpr (i == 0) return v.x;
  //   else if constexpr (i == 1) return v.y;
  //   else if constexpr (i == 2) return v.z;
  // }
  //
  // template <size_t i, typename T, int n>
  // T&& get(dpl::vector_n<T, n>&& v) {
  //        if constexpr (i == 0) return v.x;
  //   else if constexpr (i == 1) return v.y;
  //   else if constexpr (i == 2) return v.z;
  // }
}
// NOLINTEND(cert-dcl58-cpp)

// namespace std
// {
//   template <size_t i, typename T, int n>
//   decltype(auto) get(const dpl::vector_n<T, n>& v) {
//          if constexpr (i == 0) return v.x;
//     else if constexpr (i == 1) return v.y;
//     else if constexpr (i == 2) return v.z;
//   }
//
//   template <size_t i, typename T, int n>
//   decltype(auto) get(dpl::vector_n<T, n>& v) {
//          if constexpr (i == 0) return v.x;
//     else if constexpr (i == 1) return v.y;
//     else if constexpr (i == 2) return v.z;
//   }
//
//   template <size_t i, typename T, int n>
//   decltype(auto) get(dpl::vector_n<T, n>&& v) {
//          if constexpr (i == 0) return v.x;
//     else if constexpr (i == 1) return v.y;
//     else if constexpr (i == 2) return v.z;
//   }
// }


namespace dpl
{
  template <int count, int dim>
  struct sdim;

  template <int count, int dim>
  struct _sdim_impl : std::integral_constant<int, dim>
  {
    static constexpr auto next() {
      return std::conditional_t<
        dim < count - 1,
        sdim<count, dim + 1>,
        sdim<count, 0>>{};
    }
  };

  template <int count, int dim>
  struct sdim : _sdim_impl<count, dim> {};


  /**
   * \tparam face \n
   *   0 - x min\n
   *   1 - x max\n
   *   2 - y min\n
   *   3 - y max\n
   *   4 - z min\n
   *   5 - z max\n
   */
  template <int face>
  struct sface_t : std::integral_constant<int, face>
  {
    static constexpr auto dim = sdim<3, face/2>{};
    
    static constexpr auto is_upper = std::integral_constant<bool, face%2>{};

    /* inner direction */
    static constexpr auto non_zero_component = std::integral_constant<int, is_upper ? -1 : 1>{};
    static constexpr auto normal = vector3i{dim, non_zero_component};
  };

  template<int dim>
  struct cdims
  {
    cdims() = default;

    template<int face>
    explicit constexpr cdims(sface_t<face>) {}

    static constexpr auto e0 = dpl::sdim<3, dim>{};
    static constexpr auto e1 = e0.next();
    static constexpr auto e2 = e1.next();

    template <typename Tuple>
    static constexpr auto tie(Tuple& t) {
      return std::tie(t[e0], t[e1], t[e2]);
    }

    static constexpr auto i = e0;
    static constexpr auto j = e1;
    static constexpr auto k = e2;
  };

  template <int face>
  cdims(sface_t<face>) -> cdims<sface_t<face>::dim>;


  template <int dim0, int dim1> requires (dim0 != dim1 && dim0 < 3 && dim1 < 3)
  constexpr auto third(std::integral_constant<int, dim0> = {}, std::integral_constant<int, dim1> = {}) {
    return std::integral_constant<int, non_zero_idx(
      cross(
        vector3i{std::integral_constant<int, dim0>{}},
        vector3i{std::integral_constant<int, dim1>{}}
      )
    )>{};
  }


  template <typename R>
  class idx1d_map
  {
    R x_, xy_;

  public:
    idx1d_map() = default;

    template <typename T>
    idx1d_map(const vector_n<T, 3>& dim)
      : x_(dim.x), xy_(static_cast<R>(dim.x)*dim.y) {}

    R operator()(auto x, auto y, auto z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    template <typename T>
    R operator()(const vector_n<T, 3>& v) const {
      return static_cast<R>(v.x) + x_*v.y + xy_*v.z;
    }

    constexpr auto operator()(std::integral_constant<int, 0>) const { return R{1}; }
    auto operator()(std::integral_constant<int, 1>) const { return x_; }
    auto operator()(std::integral_constant<int, 2>) const { return xy_; }
  };

  template <typename T>
  idx1d_map(const vector_n<T, 3>&) -> idx1d_map<T>;
}