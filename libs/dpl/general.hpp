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

#include <cmath>
#include <fstream>
#include <future>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <vector>
#include <concepts>
#include <span>

#define def_attrib(name) \
  inline constexpr struct name##_t { \
    auto& operator()(auto& s, const auto i) const { return s(*this, i); } \
    auto& operator()(auto* s, const auto i) const { return (*s)(*this, i); } \
  } name;  // NOLINT(bugprone-macro-parentheses)


namespace dpl
{
  template <typename T>
  struct full_range_iterator
  {
    int value;

    constexpr T operator*() const noexcept {
      return T(value);
    }

    constexpr auto& operator++() noexcept {
      ++value;
      return *this;
    }

    constexpr bool operator!=(full_range_iterator rhs) const noexcept {
      return value != rhs.value;
    }
  };

  template <typename T>
  struct full_range_unsigned
  {
    using iter = full_range_iterator<T>;

    static constexpr iter begin() noexcept { return iter{0}; }
    static constexpr iter end()   noexcept { return iter{1 << sizeof(T)*8}; }
  };

  /*
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   */

  class stream_reader
  {
    std::ifstream stream_;

  public:
    explicit stream_reader(const auto& path, std::ios_base::openmode mode = std::ios::binary)
      : stream_{path, mode} {
    }

    template <typename T>
    auto& operator()(T& val) {
      stream_.read(reinterpret_cast<char*>(&val), sizeof(T));
      return *this;
    }

    template <typename T>
    auto& operator()(T* ptr, auto size) {
      stream_.read(reinterpret_cast<char*>(ptr), size*sizeof(T));
      return *this;
    }

    template <typename T>
    auto& operator()(std::unique_ptr<T[]>& ptr, auto size) {
      return (*this)(ptr.get(), size);
    }
  };

  class stream_writer
  {
    std::ofstream stream_;

  public:
    explicit stream_writer(const auto& path, std::ios_base::openmode mode = std::ios::binary)
      : stream_{path, mode} {
    }

    template <typename T>
    auto& operator()(T val) {
      stream_.write(reinterpret_cast<const char*>(&val), sizeof(T));
      return *this;
    }

    template <typename T>
    auto& operator()(T* ptr, auto size) {
      stream_.write(reinterpret_cast<const char*>(ptr), size*sizeof(T));
      return *this;
    }

    template <typename T>
    auto& operator()(const std::unique_ptr<T[]>& ptr, auto size) {
      return (*this)(ptr.get(), size);
    }
  };

  template<typename T>
  constexpr T pow(T value, std::integral auto exponent) {
    return exponent == 0 ? 1 : value*pow(value, exponent - 1);
  }

  /*
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   */

  class default_map
  {
    struct unity_multiplier
    {
      template <class T>
      constexpr T&& operator*(T&& rhs) const noexcept {
        return std::forward<T>(rhs);
      }

      template <class T>
      friend T&& operator*(T&& lhs, unity_multiplier) {
        return std::forward<T>(lhs);
      }

      // template <class T>
      // constexpr auto operator*(const T& right) const noexcept {
      //     return right;
      // }

      // template <class T>
      // constexpr const T& operator*(const T& right) const noexcept {
      //   return std::forward<T>(right);
      // }    
    };

    static inline auto true_  = [](auto) { return true; };
    static inline auto unity_ = [](auto) { return unity_multiplier{}; };
    static inline auto nan_   = [](auto) { return std::numeric_limits<double>::quiet_NaN(); };

  public:
    using true_t = decltype(true_);
    using unity_t = decltype(unity_);
    using nan_t = decltype(nan_);

    static bool invert(std::true_type, bool v) { return !v; }
    static bool invert(std::false_type, bool v) { return v; }
  };

  template <std::integral T, typename Tag, bool Validable = false, T InvalidValue = std::numeric_limits<T>::max()>
  struct strong_integer
  {
    using type = T;

    T value;

    strong_integer() noexcept = default;
    constexpr explicit strong_integer(T v) noexcept : value{v} {}

    constexpr auto& operator++() noexcept {
      ++value;
      return *this;
    }

    auto operator++(int) noexcept {
      auto temp = *this;
      ++value;
      return temp;
    }

    constexpr explicit operator T() const noexcept { return value; }

    constexpr auto& operator*() const noexcept { return value; }
              auto& operator*()       noexcept { return value; }

    constexpr bool operator< (std::integral auto rhs) const noexcept { return value <  rhs; }  // NOLINT(clang-diagnostic-sign-compare)
    constexpr bool operator==(std::integral auto rhs) const noexcept { return value == rhs; }
    constexpr auto operator+ (std::integral auto rhs) const noexcept { return strong_integer{value + rhs}; }
    constexpr auto operator- (std::integral auto rhs) const noexcept { return strong_integer{value - rhs}; }
    constexpr auto operator* (std::integral auto rhs) const noexcept { return strong_integer{value * rhs}; }

    constexpr bool operator< (strong_integer rhs) const noexcept { return value <  *rhs; }
    constexpr bool operator==(strong_integer rhs) const noexcept { return value == *rhs; }
    constexpr auto operator+ (strong_integer rhs) const noexcept { return strong_integer{value + *rhs}; }
    constexpr auto operator- (strong_integer rhs) const noexcept { return strong_integer{value - *rhs}; }

    explicit constexpr operator bool() const noexcept requires (Validable) {
      return value != InvalidValue;
    }

    static constexpr auto invalid_value() noexcept requires (Validable) {
      return strong_integer{InvalidValue};
    }
  };



  template <typename...> // template <typename Index, typename Value>
  class strong_vector {};

  template <typename U, typename Tag, bool V, U W, typename Value>
  class strong_vector<strong_integer<U, Tag, V, W>, Value>
  {
  protected:
    std::unique_ptr<Value[]> uptr_;

  public:
    strong_vector() = default;
    strong_vector(const strong_vector& other) = delete;
    strong_vector(strong_vector&& other) noexcept = default;
    strong_vector& operator=(const strong_vector& other) = delete;
    strong_vector& operator=(strong_vector&& other) noexcept = default;
    ~strong_vector() = default;

    template <std::size_t value>
    explicit constexpr strong_vector(std::integral_constant<std::size_t, value>)
      : uptr_{std::make_unique<Value[]>(value)} {}

    explicit strong_vector(strong_integer<U, Tag, V, W> size)
      : uptr_{std::make_unique<Value[]>(*size)} {}

    explicit strong_vector(strong_integer<U, Tag, V, W> size, Value value)
      : strong_vector{size} {
      std::fill_n(uptr_.get(), *size, value);
    }

    void clear() {
      uptr_.reset();
    }

    void resize(strong_integer<U, Tag, V, W> size) {
      uptr_ = std::make_unique<Value[]>(*size);
    }

    void assign(strong_integer<U, Tag, V, W> size, Value value) {
      resize(size);
      std::fill_n(uptr_.get(), *size, value);
    }

    auto& operator[](strong_integer<U, Tag, V, W> index) {
      return uptr_[*index];
    }

    auto& operator[](strong_integer<U, Tag, V, W> index) const {
      return uptr_[*index];
    }

    auto data() const {
      return uptr_.get();
    }

    explicit operator bool() const {
      return static_cast<bool>(uptr_);
    }

    auto span(const strong_integer<U, Tag, V, W> size) const {
      return std::span{data(), data() + *size};
    }
  };

  template <typename U, typename Tag, bool V, U W>
  class strong_vector<strong_integer<U, Tag, V, W>, bool>
  {
  protected:
    std::vector<bool> vec_;

  public:
    strong_vector() = default;
    strong_vector(const strong_vector& other) = delete;
    strong_vector(strong_vector&& other) noexcept = default;
    strong_vector& operator=(const strong_vector& other) = delete;
    strong_vector& operator=(strong_vector&& other) noexcept = default;
    ~strong_vector() = default;

    template <std::size_t value>
    explicit constexpr strong_vector(std::integral_constant<std::size_t, value>)
      : vec_(value) {}

    explicit strong_vector(strong_integer<U, Tag, V, W> size)
      : vec_(*size) {}

    auto operator[](strong_integer<U, Tag, V, W> i) {
      return vec_[*i];
    }

    auto operator[](strong_integer<U, Tag, V, W> i) const {
      return vec_[*i];
    }
  };


  template <typename Index, typename Value, std::size_t size = std::size_t{1} << 8*sizeof(Index)>
  class strong_array {};



  template <typename U, typename Tag, bool V, U W, typename Value, std::size_t size>
  class strong_array<strong_integer<U, Tag, V, W>, Value, size> : public strong_vector<strong_integer<U, Tag, V, W>, Value>
  {
  public:
    constexpr strong_array()
      : strong_vector<strong_integer<U, Tag, V, W>, Value>(std::integral_constant<size_t, size>{}) {}

    explicit constexpr strong_array(Value v)
      : strong_array() {
      std::fill_n(this->uptr_.get(), size, v);
    }
    
    auto begin() const {
      return this->data();
    }
    
    auto end() const {
      return this->data() + size;
    }

    // auto span() const {
    //   return std::span{this->data(), this->data() + size};
    // }
  };

  template <typename U, typename Tag, bool V, U W, std::size_t size>
  class strong_array<strong_integer<U, Tag, V, W>, bool, size> : public strong_vector<strong_integer<U, Tag, V, W>, bool>
  {
  public:
    constexpr strong_array()
      : strong_vector<strong_integer<U, Tag, V, W>, bool>(std::integral_constant<size_t, size>{}) {}
    
    auto begin() const {
      return this->vec_.begin();
    }
    
    auto end() const {
      return this->vec_.end();
    }

    // auto indices() const {
    //   return std::span<strong_integer<U, Tag, V, W>, size>{strong_integer<U, Tag, V, W>{0}, size};
    // }
  };



  template <typename... T>
  constexpr bool always_false = false;
  
  template <typename... Args>
  struct has_any : std::bool_constant<sizeof...(Args) != 0> {};

  template <typename... Args>
  inline constexpr bool has_any_v = has_any<Args...>::value;


  // template <typename V, typename... Rest>
  // struct are_arithmetic : std::bool_constant<std::is_arithmetic_v<V> && are_arithmetic<Rest...>::value> {};

  // template <typename V_>
  // struct are_arithmetic<V_> : std::is_arithmetic<V_> {};


  template <typename U, typename V, typename... Rest>
  struct are_assignable : std::bool_constant<std::is_assignable_v<U, V> && are_assignable<U, Rest...>::value> {};

  template <typename U, typename V>
  struct are_assignable<U, V> : std::is_assignable<U, V> {};

  template <typename... Args>
  inline constexpr bool are_assignable_v = are_assignable<Args...>::value;


  // template <int n, int i = 0> requires (n > 0)
  // constexpr auto tail_aggregate(const auto& op, const auto& data) {
  //   if constexpr (i == n - 1)
  //     return data[i];
  //   else
  //     return op(data[i], tail_aggregate<n, i + 1>(op, data));
  // }

  template <int n, typename T = void, int i = 0> requires (n > 0)
  constexpr auto tail_aggregate(const auto& op, const auto& data) {
    if constexpr (i == n - 1) {
      if constexpr (std::is_same_v<T, void>)
        return data[i];
      else
        return static_cast<T>(data[i]);
    }
    else
      return op(data[i], tail_aggregate<n, T, i + 1>(op, data));
  }
  
  template <int n0, int n1, typename Func>
  constexpr void sfor(const Func& f) {
    if constexpr (n0 < n1) {
      if constexpr (std::is_invocable_v<Func, std::integral_constant<int, n0>>)
        f(std::integral_constant<int, n0>{});
      else if constexpr (std::is_invocable_v<Func>)
        f();
      else
        f.template operator()<n0>();

      dpl::sfor<n0 + 1, n1>(f);
    }
  }

  template <int n, typename Func>
  constexpr void sfor(const Func& f) {
    dpl::sfor<0, n>(f);
  }

  namespace helper
  {
    template <int n0, int n1, typename Func>
    constexpr void psfor(std::future<void>* queue, const Func& f) {
      if constexpr (n0 < n1) {
        if constexpr (std::is_invocable_v<Func, std::integral_constant<int, n0>>)
          *queue = std::async(std::launch::async, f, std::integral_constant<int, n0>{});
        else if constexpr (std::is_invocable_v<Func>)
          *queue = std::async(std::launch::async, f);
        // else
        //   *queue = std::async(std::launch::async, [&f] { f.template operator()<n0>(); });

        helper::psfor<n0 + 1, n1>(++queue, f);
      }
    }
  }
  
  template <int n, typename Func>
  constexpr void psfor(const Func& f) {
    std::array<std::future<void>, n> queue;
    
    helper::psfor<0, n>(queue.data(), f);

    sfor<n>([&](auto i) {
      queue[i].wait();
    });
    
    // for (std::future<void>& task : queue)
    //   task.wait();
  }

  namespace extrapolant
  {
    static inline constexpr struct linear_t {} linear;
    static inline constexpr struct flat_t {} flat;
  }

  namespace lerp_base
  {
    static inline constexpr struct linear_t {} linear;
    static inline constexpr struct log10_t {} log10;
  };
  
  template <typename Base = lerp_base::linear_t, typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(const T& a, const T& b, double t) {
    if constexpr (std::is_same_v<Extrapolant, extrapolant::linear_t> && std::is_same_v<Base, lerp_base::linear_t>)
      return a + (b - a)*t;
    else if constexpr (std::is_same_v<Base, lerp_base::log10_t>)
      return std::pow(10., dpl::lerp<lerp_base::linear_t, Extrapolant>(std::log10(a), std::log10(b), t));
    else if constexpr (std::is_same_v<Extrapolant, extrapolant::flat_t>)
      return
        t < 0 ? a :
        t > 1 ? b :
        dpl::lerp<Base, extrapolant::linear_t>(a, b, t);
    else
      static_assert(always_false<T>, "wrong extrapolant or base");
  }

  template <typename Base = lerp_base::linear_t, typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(T* ptr, double t) {
    return dpl::lerp<Base, Extrapolant>(ptr[0], ptr[1], t);
  }

  template <typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(bool log10, const T& a, const T& b, double t) {
    return log10
      ? dpl::lerp<lerp_base::log10_t, Extrapolant>(a, b, t)
      : dpl::lerp<lerp_base::linear_t, Extrapolant>(a, b, t);
  }

  template <typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(bool log10, T* ptr, double t) {
    return dpl::lerp<Extrapolant>(log10, ptr[0], ptr[1], t);
  }
  


  template <char ThousandsSep = ',', typename T>
  std::string ft_format(const T& value, int precision = 0) {
    struct separate_thousands : std::numpunct<char> {
      char_type do_thousands_sep() const override { return ThousandsSep; }  // separate with commas
      string_type do_grouping() const override { return "\3"; } // groups of 3 digit
    };

    std::ostringstream ss;
    ss.imbue(std::locale(ss.getloc(), new separate_thousands));
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
  }
  

 
}


template <std::integral U, typename Tag, bool V, U W>
struct std::incrementable_traits<dpl::strong_integer<U, Tag, V, W>> : std::incrementable_traits<U> {};  // NOLINT(cert-dcl58-cpp)

template <std::integral U, typename Tag, bool V, U W>
struct std::hash<dpl::strong_integer<U, Tag, V, W>>  // NOLINT(cert-dcl58-cpp)
{
  std::size_t operator()(const dpl::strong_integer<U, Tag, V, W>& k) const {
    return std::hash<U>()(*k);
  }
};

// namespace std
// {
//   // template <integral T, typename Tag>
//   // struct incrementable_traits<dpl::strong_integer<U, Tag, V, W>> : incrementable_traits<T> {};
//
//   // template <integral T, typename Tag>
//   // struct hash<dpl::strong_integer<U, Tag, V, W>>
//   // {
//   //   size_t operator()(const dpl::strong_integer<U, Tag, V, W>& k) const {
//   //     return hash<T>()(*k);
//   //   }
//   // };
// }
