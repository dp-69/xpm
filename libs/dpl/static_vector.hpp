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

#include <ostream>
#include <span>

namespace dpl
{
  static inline constexpr auto _0 = std::integral_constant<int, 0>{};

  struct operation
  {
    static constexpr auto add = [](const auto& l, const auto& r) { return l + r; };
    static constexpr auto diff = [](const auto& l, const auto& r) { return l - r; };
    static constexpr auto mult = [](const auto& l, const auto& r) { return l*r; };
    static constexpr auto div = [](const auto& l, const auto& r) { return l/r; };
    static constexpr auto min = [](const auto& l, const auto& r) { return std::min(l, r); };
    
    struct cross
    {
      static constexpr auto x = [](const auto& a, const auto& b) {
        return a[1]*b[2] - a[2]*b[1]; 
      };
      
      static constexpr auto y = [](const auto& a, const auto& b) {
        return a[2]*b[0] - a[0]*b[2]; 
      };

      static constexpr auto z = [](const auto& a, const auto& b) {
        return a[0]*b[1] - a[1]*b[0]; 
      };
    };
  };
  
  
  template <typename T_, int n_>
  class vector_n;
  

  template <typename Derived_, typename T_, int i_>
  class _vector_rec : public _vector_rec<Derived_, T_, i_ - 1> {};

  template <typename Derived_, typename T_>
  class _vector_rec<Derived_, T_, 0>
  {
  public:
    auto& x() { return static_cast<Derived_*>(this)->template _get<0>(); }
    constexpr auto& x() const { return static_cast<const Derived_*>(this)->template _get<0>(); }
  };

  template <typename Derived_, typename T_>
  class _vector_rec<Derived_, T_, 1> : public _vector_rec<Derived_, T_, 0>
  {
  public:
    auto& y() { return static_cast<Derived_*>(this)->template _get<1>(); }
    constexpr auto& y() const { return static_cast<const Derived_*>(this)->template _get<1>(); }

    constexpr auto slope() const { return this->y()/this->x(); }
  };
  
  template <typename Derived_, typename T_>
  class _vector_rec<Derived_, T_, 2> : public _vector_rec<Derived_, T_, 1>
  {        
  public:
    auto& z() { return static_cast<Derived_*>(this)->template _get<2>(); }
    constexpr auto& z() const { return static_cast<const Derived_*>(this)->template _get<2>(); }

    // template <typename U_ = vector_n<T_, 3>>
    // constexpr auto cross(const U_& r) const {
    //   return vector_n<T_, 3>{operation::cross{}, *static_cast<const Derived_*>(this), r};
    // }
  };


  template <typename Derived_, typename T_, int n_>
  class _vector_oper : public _vector_rec<Derived_, T_, n_ - 1>
  {
    template <typename U_, typename Func_>
    constexpr auto _binary_op_impl_ptr(U_* r, const Func_& f) const {
      vector_n<decltype(f(T_{}, U_{})), n_> v;
      sfor<n_>([&](auto i) { v[i] = f((*this)[i], r[i]); });
      return v;
    }

    // template <typename U_, typename Func_>
    // constexpr auto _binary_op_impl_const(U_ r, const Func_& f) const {
    //   // vector_n<decltype(f(T_{}, U_{})), n_> v;
    //   // s_for<n_>([&](auto i) { v[i] = f((*this)[i], r); });
    //   // return v;
    //   return vector_n<decltype(f(T_{}, U_{})), n_>{*this, r, f};
    // }

    // template <typename U_, typename Func_>
    // constexpr auto _binary_op_impl_vec(const U_& r, const Func_& f) const {
    //   // vector_n<decltype(f(T_{}, r[_0])), n_> v{};
    //   // s_for<n_>([&](auto i) { v[i] = f((*this)[i], r[i]); });
    //   // return v;
    //   return vector_n<decltype(f(T_{}, r[_0])), n_>{*this, r, f};
    // }

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
    
    template <int i_>
    auto& get() {
      return static_cast<Derived_*>(this)->template _get<i_>();
    }

    template <int i_>
    constexpr auto& get() const {
      return static_cast<const Derived_*>(this)->template _get<i_>();
    }


    // explicit operator std::string() {
    //   std::stringstream ss;
    //   ss << *this;
    //   return ss.str();
    // }
    
    template <typename U_, U_ i_>
    auto& operator[](std::integral_constant<U_, i_>) {
      // static_assert(false, "ERROR");
      return get<i_>();
    }

    template <typename U_, U_ i_>
    constexpr auto& operator[](std::integral_constant<U_, i_>) const {
      // static_assert(false, "ERROR");
      return get<i_>();
    }

    template <typename U_> requires std::is_integral_v<U_>
    auto& operator[](U_ i) {
      // static_assert(false, "ERROR");
      return static_cast<Derived_*>(this)->_get(i);
    }
    
    template <typename U_> requires std::is_integral_v<U_>
    constexpr auto& operator[](U_ i) const {
      // static_assert(false, "ERROR");
      return static_cast<const Derived_*>(this)->_get(i);
    }


    template <typename U_ = T_>
    constexpr auto prod() const {
      U_ val = (*this)[_0];
      sfor<1, n_>([&](auto i) { val *= (*this)[i]; });
      return val;
    }

    template <typename U_ = T_>
    constexpr auto sum() const {
      U_ val = (*this)[_0];
      sfor<1, n_>([&](auto i) { val += (*this)[i]; });
      return val;
    }

    // template <typename U_ = T_>
    // constexpr auto min() const {
    //   U_ val = (*this)[_0];
    //   s_for<1, n_>([&](auto i) { val = std::min(val, (*this)[i]); });
    //   return val;
    // }

    template <typename U_ = vector_n<T_, n_>>
    constexpr auto dot(const U_& r) const {
      return (*this*r).sum();
      // _binary_op_impl(r, mult_).sum();
    }

    // template <typename U_>
    // constexpr auto cross(const U_& r) const {
    //   return vector_n<T_, n_>{operation::cross{}, *this, r};
    // }

    constexpr auto length() const {
      return std::sqrt(this->dot(*this));
    }
    
    constexpr auto normalise() const {
      return *this/this->length();
    }
    
    template <typename U_>
    constexpr auto operator+(const U_& r) const { return _binary_op_impl(r, operation::add); }

    template <typename U_>
    constexpr auto operator-(const U_& r) const { return _binary_op_impl(r, operation::diff); }

    template <typename U_>
    constexpr auto operator*(const U_& r) const { return _binary_op_impl(r, operation::mult); }

    template <typename U_>
    constexpr auto operator/(const U_& r) const { return _binary_op_impl(r, operation::div); }

    constexpr auto min() const {
      auto min = (*this)[0];
      
      for (auto i = 1; i < n_; ++i)
        min = std::min(min, (*this)[i]);

      return min;
    }

    template <typename U_ = vector_n<T_, n_>>
    constexpr auto min(const U_& r) const {
      return _binary_op_impl(r, operation::min);
    }

    constexpr auto max() const {
      auto min = (*this)[0];
      
      for (auto i = 1; i < n_; ++i)
        min = std::max(min, (*this)[i]);

      return min;
    }
  };

  template <typename T, int i>
  constexpr auto cross(const vector_n<T, i>& l, const vector_n<T, i>& r) {
    return vector_n<T, i>{operation::cross{}, l, r};
  }

  template <typename T, int n>
  constexpr int non_zero_idx(const vector_n<T, n>& v) {
    for (auto i = 0; i < n; ++i)
      if (v[i])
        return i;

    return -1;
  }


  template <typename Derived, int n>
  constexpr auto operator/(double l, const _vector_oper<Derived, double, n>& r) {
    vector_n<double, n> v;
    sfor<n>([&](auto i) { v[i] = l/r[i]; });
    return v;
  }
  
  template <typename Derived_, typename T_, int n_>
  std::ostream& operator<<(std::ostream& os, const _vector_oper<Derived_, T_, n_>& v) {
    os << v[0];
    sfor<1, n_>([&](auto i) { os << ", " << v[i]; });
    return os;
  }
  
  template <typename T, int n>
  class vector_n : public _vector_oper<vector_n<T, n>, T, n>
  {
    T ptr_[n];

    template <typename V, typename... Rest>
    constexpr void _assign(const V& val, Rest... args) {
      _get<n - sizeof...(Rest) - 1>() = val;
      if constexpr (has_any_v<Rest...>)
        _assign<Rest...>(args...);
    }

    template <int i>
    constexpr auto& _get() {
      return ptr_[i];
    }

    template <int i>
    constexpr auto& _get() const {
      return ptr_[i];
    }

    constexpr auto& _get(int i) {
      return ptr_[i];
    }

    constexpr auto& _get(int i) const {
      return ptr_[i];
    }

    friend class _vector_oper<vector_n, T, n>;
    friend class _vector_rec<vector_n, T, 0>;
    friend class _vector_rec<vector_n, T, 1>;
    friend class _vector_rec<vector_n, T, 2>;

  public:
    operator T*() { return ptr_; }
    operator const T*() const { return ptr_; }
    
    vector_n() = default;

    template <typename U> requires (std::is_same_v<U, T> || std::is_arithmetic_v<U>)
    constexpr vector_n(const U val) {
      sfor<n>([&](auto i) { ptr_[i] = val; });
    }

    template <typename ...Args> requires (sizeof...(Args) == n && are_assignable_v<T&, Args...>)
    constexpr vector_n(Args... args) {
      _assign<Args...>(args...);
    }

    // template <typename Func> requires std::is_invocable_v<Func, std::integral_constant<int, 0>>
    // explicit constexpr vector_n(const Func& f)/* : ptr_{}*/ {
    //   sfor<n>([&](auto i) { ptr_[i] = f(i); });
    // }

    template <typename U, typename V, typename Func>
      requires std::is_invocable_v<Func, typename U::ScalarType, typename V::ScalarType>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, n>& l,
      const _vector_oper<V, typename V::ScalarType, n>& r) {
      sfor<n>([&](auto i) { ptr_[i] = f(l[i], r[i]); });
    }

    template <typename U, typename V, typename Func> requires std::is_arithmetic_v<V>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, n>& l,
      const V r) {
      sfor<n>([&](auto i) { ptr_[i] = f(l[i], r); });
    }

    template <typename Derived, typename U>
    constexpr vector_n(const _vector_oper<Derived, U, n>& v) {
      sfor<n>([&](auto i) { ptr_[i] = static_cast<T>(v[i]); });
    }
    
    template <typename U, typename V> requires (n == 3)
    constexpr vector_n(operation::cross, const U& l, const V& r) {
      ptr_[0] = operation::cross::x(l, r);
      ptr_[1] = operation::cross::y(l, r);
      ptr_[2] = operation::cross::z(l, r);
    }

    template <int i, typename U = int>
    explicit constexpr vector_n(std::integral_constant<int, i>, const U v = 1) {
      sfor<n>([v, this]<int j> {
        if constexpr (i == j)
          ptr_[j] = v;
        else
          ptr_[j] = 0;
      });
    }
    
    vector_n(auto* ptr) {
      sfor<n>([&](auto i) { ptr_[i] = ptr[i]; });
    }

    vector_n(const auto* ptr) {
      sfor<n>([&](auto i) { ptr_[i] = ptr[i]; });
    }
  };






  
  template <typename T, int n>
  class vector_n_map : public _vector_oper<vector_n_map<T, n>, T, n>
  {
    T* ptr_;

    template <int i>
    auto& _get() {
      return ptr_[i];
    }

    template <int i>
    const auto& _get() const {
      return ptr_[i];
    }

    auto& _get(int i) {
      return ptr_[i];
    }

    const auto& _get(int i) const {
      return ptr_[i];
    }

    friend class _vector_oper<vector_n_map, T, n>;
    friend class _vector_rec<vector_n_map, T, 0>;
    friend class _vector_rec<vector_n_map, T, 1>;
    friend class _vector_rec<vector_n_map, T, 2>;

  public:
    operator T*() { return ptr_; }
    operator const T*() const { return ptr_; }

    vector_n_map(T* ptr) : ptr_(ptr) {}
  };
   


  using vector2i = vector_n<int, 2>;
  using vector2d = vector_n<double, 2>;
  using vector2i_map = vector_n_map<int, 2>;
  using vector3i = vector_n<int, 3>;
  using vector3d = vector_n<double, 3>;
  using vector3i_map = vector_n_map<int, 3>;
  using vector3d_map = vector_n_map<double, 3>;



  template <typename T>
  auto solve(const vector_n<T, 2>& p0, const vector_n<T, 2>& p1, const auto arg) {
    return p0.y() + (p1 - p0).slope()*(arg - p0.x());
  }

  template <typename T>
  auto solve(std::span<const vector_n<T, 2>> curve, const auto arg, extrapolant::linear_t = extrapolant::linear) {
    if (arg <= curve.front().x())
      return solve(curve[0], curve[1], arg);

    if (arg > curve.back().x())
      return solve(curve[curve.size() - 2], curve[curve.size() - 1], arg);

    auto iter = std::ranges::lower_bound(curve, arg, {}, [](const vector_n<T, 2>& p) { return p.x(); });
    return solve(*(iter - 1), *iter, arg);
  }

  template <typename T>
  auto solve(std::span<const vector_n<T, 2>> curve, const auto arg, extrapolant::flat_t) {
    if (arg <= curve.front().x())
      return curve.front().y();

    if (arg > curve.back().x())
      return curve.back().y();

    auto iter = std::ranges::lower_bound(curve, arg, {}, [](const vector_n<T, 2>& p) { return p.x(); });
    return solve(*(iter - 1), *iter, arg);
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

  template <size_t i, typename Type, int n>
  struct tuple_element<i, dpl::vector_n<Type, n>> {
    using type = Type;
  };

  template <typename Type, int n>
  struct tuple_size<dpl::vector_n_map<Type, n>> : integral_constant<size_t, n> {};    

  template <size_t i, typename Type, int n>
  struct tuple_element<i, dpl::vector_n_map<Type, n>> {
    using type = Type;
  };
}
// NOLINTEND(cert-dcl58-cpp)


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
  struct sface : std::integral_constant<int, face>
  {
    static constexpr auto dim = sdim<3, face/2>{};
    
    static constexpr auto is_upper = std::integral_constant<bool, face%2>{};

    /* inner direction */
    static constexpr auto non_zero_component = std::integral_constant<int, is_upper ? -1 : 1>{};
    static constexpr auto normal = vector3i{dim, non_zero_component};
    /* */
  };

  template<int dim>
  struct cdims
  {
    cdims() = default;

    template<int face>
    explicit constexpr cdims(sface<face>) {}

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
  cdims(sface<face>) -> cdims<sface<face>::dim>;


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
      : x_(dim.x()), xy_(static_cast<R>(dim.x())*dim.y()) {}

    R operator()(auto x, auto y, auto z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    template <typename T>
    R operator()(const vector_n<T, 3>& v) const {
      return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    }

    auto operator()(std::integral_constant<int, 0>) const { return 1; }
    auto operator()(std::integral_constant<int, 1>) const { return x_; }
    auto operator()(std::integral_constant<int, 2>) const { return xy_; }
  };

  template <typename T>
  idx1d_map(const vector_n<T, 3>&) -> idx1d_map<T>;
}