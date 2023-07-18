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

#include <dpl/general.hpp>

#include <ostream>

// TODO includes
#include <nlohmann/json.hpp>
// #include <QColor>


// #include <tuple>

namespace dpl
{
  static inline constexpr auto _0 = ic<0>{};
  static inline constexpr auto _1 = ic<1>{};
  static inline constexpr auto _2 = ic<2>{};

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

    template <typename U_ = vector_n<T_, 3>>
    constexpr auto cross(const U_& r) const {
      return vector_n<T_, 3>{operation::cross{}, *static_cast<const Derived_*>(this), r};
    }
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

    template <typename U_, std::enable_if_t<std::is_integral_v<U_>, int> = 0>
    auto& operator[](U_ i) {
      // static_assert(false, "ERROR");
      return static_cast<Derived_*>(this)->_get(i);
    }
    
    template <typename U_, std::enable_if_t<std::is_integral_v<U_>, int> = 0>
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

    constexpr int non_zero_dim() const {
      for (auto i = 0; i < n_; ++i)
        if ((*this)[i])
          return i;

      return -1;
    }

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

    
    // template<int _nn = n_, std::enable_if_t<_nn == 2 && std::is_floating_point_v<T_>, int> = 0>
    // auto slope() {
    //   return this->y()/this->x();
    // }
    
    // static constexpr bool is_two_element = n_ == 2;
    
    // template </*typename U_ = T_, */std::enable_if_t<true
    // /* && std::is_floating_point_v<U_>*/
    // /*std::is_same_v<n_, 2> && std::is_floating_point_v<T_>*/, int> = 0>
    // double
    // std::enable_if_t<std::is_same_v<n_, 2> && std::is_floating_point_v<T_>, T_>
    // std::enable_if_t<is_two_element, double>
    // slope()  {
      // return 3;
      // return this->y()/this->x();    
    // }



    // template <typename U_, std::enable_if_t<std::is_floating_point_v<U_>, int> = 0>
    // operator QColor() const {
    //   Qt::GlobalColor::white
    //   // static_assert(false, "ERROR");
    //   return static_cast<const Derived_*>(this)->_get(i);
    // }

    
  };



  template<typename Derived, int n>
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



  
#ifdef INCLUDE_NLOHMANN_JSON_HPP_
  template <typename Derived, typename T, int n>
  void to_json(nlohmann::json& j, const _vector_oper<Derived, T, n>& v) {
    sfor<n>([&](auto i) { j.push_back(v[i]); });
  }

  template <typename Derived, typename T, int n>
  void from_json(const nlohmann::json& j, _vector_oper<Derived, T, n>& v) {
    sfor<n>([&](auto i) { v[i] = j[i]; });
  }
#endif
  






  
  template <typename T, int n>
  class vector_n : public _vector_oper<vector_n<T, n>, T, n>
  {
    T ptr_[n];

    template <typename V_, typename ...Rest>
    constexpr void _assign(const V_& val, Rest... args) {
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

    friend class _vector_oper<vector_n<T, n>, T, n>;
    friend class _vector_rec<vector_n<T, n>, T, 0>;
    friend class _vector_rec<vector_n<T, n>, T, 1>;
    friend class _vector_rec<vector_n<T, n>, T, 2>;

  public:
    operator T*() { return ptr_; }
    operator const T*() const { return ptr_; }
    
    vector_n() = default;


    template <typename U, std::enable_if_t<std::is_same_v<U, T> || std::is_arithmetic_v<U>, int> = 0>
    constexpr vector_n(const U val) : ptr_{} {
      sfor<n>([&](auto i) { ptr_[i] = val; });
    }

    template <typename ...Args,
      std::enable_if_t<sizeof...(Args) == n && 
      are_assignable<T&, Args...>::value
      , int> = 0>
    constexpr vector_n(Args... args) : ptr_{} {
      // static_assert(sizeof...(Args) == n_, "Number of arguments is not equal to the size of vector");
      _assign<Args...>(args...);
    }

    template <typename U, typename V, typename Func>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, n>& l,
      const _vector_oper<V, typename V::ScalarType, n>& r) : ptr_{} {
      sfor<n>([&](auto i) { ptr_[i] = f(l[i], r[i]); });
    }

    template <typename U, typename V, typename Func,
      std::enable_if_t<std::is_arithmetic_v<V>, int> = 0>
    constexpr vector_n(
      const Func& f,
      const _vector_oper<U, typename U::ScalarType, n>& l,
      const V& r) : ptr_{} {
      sfor<n>([&](auto i) { ptr_[i] = f(l[i], r); });
    }

    template <typename Derived, typename U>
    constexpr vector_n(const _vector_oper<Derived, U, n>& v) : ptr_{} {
      sfor<n>([&](auto i) { ptr_[i] = v[i]; });
    }
    
    
    template <typename U, typename V>
    constexpr vector_n(
      operation::cross,
      const U& l,
      const V& r) : ptr_{} {

      ptr_[0] = operation::cross::x(l, r);
      ptr_[1] = operation::cross::y(l, r);
      ptr_[2] = operation::cross::z(l, r);
    }

    template <int i, typename U = int>
    constexpr vector_n(ic<i>, const U& v = 1) : ptr_{} {
      ptr_[i] = v;
    }

    

    
    template <typename U>
    vector_n(U* ptr) {
      sfor<n>([&](auto i) { ptr_[i] = ptr[i]; });
    }

    template <typename U>
    vector_n(const U* ptr) {
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

    friend class _vector_oper<vector_n_map<T, n>, T, n>;
    friend class _vector_rec<vector_n_map<T, n>, T, 0>;
    friend class _vector_rec<vector_n_map<T, n>, T, 1>;
    friend class _vector_rec<vector_n_map<T, n>, T, 2>;

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


  template<typename T, typename U>
  auto lerp_solve(const vector_n<T, 2>& p0, const vector_n<T, 2>& p1, U arg) {
    return p0.y() + (p1 - p0).slope()*(arg - p0.x());
  }
  
  template<typename T, typename U>
  auto lerp_solve(const std::vector<vector_n<T, 2>>& curve, U arg) {
    auto begin = curve.begin();
    auto end = curve.end();

    if (arg < begin->x())
      return lerp_solve(*begin, *(begin + 1), arg);        

    auto iter = std::lower_bound(begin, end, arg,
      [](const vector_n<T, 2>& l, T r) { return l.x() < r; });

    if (iter == end - 1 || iter == end)
      return lerp_solve(*(end - 2), *(end - 1), arg);
    
    return lerp_solve(*iter, *(iter + 1), arg);
  }
}

namespace std
{
  template <typename Type, int n>
  struct tuple_size<dpl::vector_n<Type, n>> : integral_constant<size_t, n> {};        

  template <size_t i, typename Type, int n>
  struct tuple_element<i, dpl::vector_n<Type, n>>
  {
    using type = Type;
  };


  template <typename Type, int n>
  struct tuple_size<dpl::vector_n_map<Type, n>> : integral_constant<size_t, n> {};    

  template <size_t i, typename Type, int n>
  struct tuple_element<i, dpl::vector_n_map<Type, n>>
  {
    using type = Type;
  };
}







namespace dpl
{
  template <int count, int dim>
  struct sdim;

  template<int count, int dim>
  struct _sdim_impl : ic<dim>
  {
    static constexpr auto next() {
      return std::conditional_t<
        dim < count - 1,
        sdim<count, dim + 1>,
        sdim<count, 0>>{};
    }
  };

  template<int count, int dim>
  struct sdim : _sdim_impl<count, dim> {};

  template<int dim>
  struct sdim<3, dim> : _sdim_impl<3, dim>
  {
    template<int dim1/*, std::enable_if_t<dim_idx1 != dim_idx, int> = 0*/>
    static constexpr auto cross(ic<dim1> = {}) {
      static_assert(dim1 != dim, "Dimensions should not be equal");
      return sdim<3, vector3i{ic<dim>{}}.cross({ic<dim1>{}}).non_zero_dim()>{};
    }
  };





  template<int dim>
  struct cdims
  {
    static constexpr auto e0 = dpl::sdim<3, dim>{};
    static constexpr auto e1 = e0.next();
    static constexpr auto e2 = e1.next();

    template<typename Tuple>
    static constexpr auto tie(Tuple& t) {
      return std::tie(t[e0], t[e1], t[e2]);
    }
  };

    

  // namespace color
  // {
  //   static inline constexpr auto white = vector3d{1.0};
  //   static inline constexpr auto black = vector3d{0.0};
  // }

  

  /**
   * \tparam face in [0, 5]
   *              x min - 0
   *              x max - 1
   *              y min - 2
   *              y max - 3
   *              z min - 4
   *              z max - 5
   */
  template<int face>
  struct face_cubic : ic<face>
  {
    static constexpr auto dim = sdim<3, face/2>{};
    
    static constexpr auto is_upper = std::integral_constant<bool, face%2>{};

    /**
     * \brief inner direction
     */
    static constexpr auto non_zero_component = std::conditional_t<is_upper, ic<-1>, ic<1>>{};

    /**
     * \brief inner direction
     */
    static constexpr auto normal = vector3i{dim, non_zero_component};
  };
}







