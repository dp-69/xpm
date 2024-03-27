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
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <cmath>
#include <future>
#include <vector>

// #ifdef __cpp_lib_ranges
//   #include <ranges>
// #endif

// #ifdef __cpp_lib_span
// #endif

#include <span>
#include <concepts>

#define def_attrib(name) \
  inline constexpr struct name##_t { \
    auto& operator()(auto& s, const auto i) const { return s(*this, i); } \
    auto& operator()(auto* s, const auto i) const { return (*s)(*this, i); } \
  } name;  // NOLINT(bugprone-macro-parentheses)


namespace dpl
{
  class default_map
  {
    static inline auto true_ = [](auto) { return true; };
    static inline auto unity_ = [](auto) { return 1; };
    static inline auto nan_ = [](auto) { return std::numeric_limits<double>::quiet_NaN(); };

  public:
    using true_t = decltype(true_);
    using unity_t = decltype(unity_);
    using nan_t = decltype(nan_);

    static bool invert(std::true_type, bool v) { return !v; }
    static bool invert(std::false_type, bool v) { return v; }
  };

  template <std::integral T, typename Tag = void>
  struct strong_integer
  {
    using type = T;
    // using difference_type = int;

    type value;

    strong_integer() = default;
    constexpr explicit strong_integer(type v) : value{v} {}
    constexpr explicit operator type() const { return value; }

    constexpr auto& operator++() {
      ++value;
      return *this;
    }

    auto operator++(int) {
      auto temp = *this;
      ++*this;
      return temp;
    }

    constexpr auto& operator*() const { return value; }
    auto& operator*() { return value; }

    constexpr bool operator<(std::integral auto rhs) const { return value < rhs; }  // NOLINT(clang-diagnostic-sign-compare)
    constexpr bool operator>=(std::integral auto rhs) const { return value >= rhs; }
    constexpr bool operator==(std::integral auto rhs) const { return value == rhs; }
    constexpr bool operator<(const strong_integer& rhs) const { return value < *rhs; }

    constexpr auto operator+(std::integral auto rhs) const { return strong_integer{value + rhs}; }
    constexpr auto operator-(std::integral auto rhs) const { return strong_integer{value - rhs}; }
    constexpr auto operator*(std::integral auto rhs) const { return strong_integer{value * rhs}; }

    template <std::integral V>
    constexpr auto operator+(const strong_integer<V, Tag>& rhs) const { return strong_integer{value + *rhs}; }

    template <std::integral V>
    constexpr auto operator-(const strong_integer<V, Tag>& rhs) const { return strong_integer{value - *rhs}; }

    constexpr bool operator==(const strong_integer& rhs) const { return value == *rhs; }
    constexpr bool operator!=(const strong_integer& rhs) const { return value != *rhs; }

    friend constexpr bool operator<(std::integral auto lhs, const strong_integer& rhs) { return lhs < *rhs; }

    // friend constexpr bool operator==(const strong_integer& lhs, const strong_integer& rhs) { return *lhs == *rhs; }
    // friend constexpr bool operator!=(const strong_integer& lhs, const strong_integer& rhs) { return *lhs != *rhs; }
  };

  


  // template <typename>
  // struct strong_traits {};
  //
  // template <std::integral T>
  // struct strong_traits<T> {
  //   static auto& get(const T& x) { return x; }
  //   static auto& get(T& x) { return x; }
  // };
  //
  // template <std::integral T, typename Tag>
  // struct strong_traits<strong_integer<T, Tag>> {
  //   static auto& get(const strong_integer<T, Tag>& x) { return *x; }
  //   static auto& get(strong_integer<T, Tag>& x) { return *x; }
  // };


  template <typename...> // template <typename Index, typename Value>
  class strong_vector {};

  template <typename T, typename Tag, typename Value>
  class strong_vector<strong_integer<T, Tag>, Value>
  {
  protected:
    std::unique_ptr<Value[]> uptr_;

  public:
    strong_vector() = default;
    strong_vector(const strong_vector& other) = delete;
    strong_vector(strong_vector&& other) noexcept = default;
    strong_vector& operator=(const strong_vector& other) = delete;
    strong_vector& operator=(strong_vector&& other) noexcept = default;

    template<std::size_t value>
    explicit constexpr strong_vector(std::integral_constant<std::size_t, value>)
      : uptr_{std::make_unique<Value[]>(value)} {}

    explicit strong_vector(strong_integer<T, Tag> size)
      : uptr_{std::make_unique<Value[]>(*size)} {}

    explicit strong_vector(strong_integer<T, Tag> size, Value value)
      : strong_vector{size} {
      std::fill_n(uptr_.get(), *size, value);
    }

    void resize(strong_integer<T, Tag> size) {
      uptr_ = std::make_unique<Value[]>(*size);
    }

    void resize(strong_integer<T, Tag> size, Value value) {
      resize(size);
      std::fill_n(uptr_.get(), *size, value);
    }

    auto& operator[](strong_integer<T, Tag> index) {
      return uptr_[*index];
    }

    auto& operator[](strong_integer<T, Tag> index) const {
      return uptr_[*index];
    }

    auto data() const {
      return uptr_.get();
    }

    explicit operator bool() const {
      return static_cast<bool>(uptr_);
    }

    auto span(strong_integer<T, Tag> size) const {
      return std::span{data(), data() + *size};
    }
  };

  // template <typename T, typename Tag, typename ValueType>
  // class strong_vector<strong_integer<T, Tag>, ValueType>
  // {
  //   // std::unique_ptr<ValueType[]> uptr_;
  //   std::vector<ValueType> vec_;
  //
  //
  // public:
  //   strong_vector() = default;
  //   strong_vector(const strong_vector& other) = delete;
  //   strong_vector(strong_vector&& other) noexcept = default;
  //   strong_vector& operator=(const strong_vector& other) = delete;
  //   strong_vector& operator=(strong_vector&& other) noexcept = default;
  //
  //   template<std::size_t value>
  //   explicit constexpr strong_vector(std::integral_constant<std::size_t, value>)
  //     : vec_(value) {}
  //     // : uptr_{std::make_unique<ValueType[]>(value)} {}
  //
  //   explicit strong_vector(strong_integer<T, Tag> size)
  //     : vec_(*size) {}
  //     // : uptr_{std::make_unique<ValueType[]>(*size)} {}
  //
  //   explicit strong_vector(strong_integer<T, Tag> size, ValueType value)
  //     // : strong_vector{size}
  //   {
  //     vec_(*size, value);
  //     // std::fill_n(uptr_.get(), *size, value);
  //   }
  //
  //   void resize(strong_integer<T, Tag> size) {
  //     vec_.resize(*size);
  //     // uptr_ = std::make_unique<ValueType[]>(*size);
  //   }
  //
  //   void resize(strong_integer<T, Tag> size, ValueType value) {
  //     vec_.resize(*size, value);
  //     // resize(size);
  //     // std::fill_n(uptr_.get(), *size, value);
  //   }
  //
  //   auto& operator[](strong_integer<T, Tag> i) {
  //     return vec_[*i];
  //     // return uptr_[*i];
  //   }
  //
  //   auto& operator[](strong_integer<T, Tag> i) const {
  //     return vec_[*i];
  //     // return uptr_[*i];
  //   }
  //
  //   auto data() const {
  //     return vec_.data();
  //   }
  //
  //   auto data() {
  //     return vec_.data();
  //   }
  //
  //   auto begin() const {
  //     return vec_.begin();
  //   }
  //
  //   auto end() const {
  //     return vec_.end();
  //   }
  // };

  template <typename T, typename Tag>
  class strong_vector<strong_integer<T, Tag>, bool>
  {
  protected:
    std::vector<bool> vec_;

  public:
    strong_vector() = default;
    strong_vector(const strong_vector& other) = delete;
    strong_vector(strong_vector&& other) noexcept = default;
    strong_vector& operator=(const strong_vector& other) = delete;
    strong_vector& operator=(strong_vector&& other) noexcept = default;

    template<std::size_t value>
    explicit constexpr strong_vector(std::integral_constant<std::size_t, value>)
      : vec_(value) {}

    explicit strong_vector(strong_integer<T, Tag> size)
      : vec_(*size) {}

    auto operator[](strong_integer<T, Tag> i) {
      return vec_[*i];
    }

    auto operator[](strong_integer<T, Tag> i) const {
      return vec_[*i];
    }
  };


  template <typename Index, typename Value, std::size_t size = std::size_t{1} << 8*sizeof(Index)>
  class strong_array {};

  template <typename T, typename Tag, typename Value, std::size_t size>
  class strong_array<strong_integer<T, Tag>, Value, size> : public strong_vector<strong_integer<T, Tag>, Value>
  {
  public:
    constexpr strong_array()
      : strong_vector<strong_integer<T, Tag>, Value>(std::integral_constant<size_t, size>{}) {}

    constexpr strong_array(Value v)
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

  template <typename T, typename Tag, std::size_t size>
  class strong_array<strong_integer<T, Tag>, bool, size> : public strong_vector<strong_integer<T, Tag>, bool>
  {
  public:
    constexpr strong_array()
      : strong_vector<strong_integer<T, Tag>, bool>(std::integral_constant<size_t, size>{}) {}
    
    auto begin() const {
      return this->vec_.begin();
    }
    
    auto end() const {
      return this->vec_.end();
    }

    // auto indices() const {
    //   return std::span<strong_integer<T, Tag>, size>{strong_integer<T, Tag>{0}, size};
    // }
  };



  template <typename... T>
  constexpr bool always_false = false;
  
  template <typename... Args_>
  struct has_any : std::bool_constant<sizeof...(Args_) != 0> {};

  template <typename... Args_>
  inline constexpr bool has_any_v = has_any<Args_...>::value;

  template <int i_>
  using ic = std::integral_constant<int, i_>;


  
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





  // std::list<std::future<void>> queue;
  


  

  template <int n0, int n1, typename Functor>
  constexpr void sfor(const Functor& f) {
    if constexpr (n0 < n1) {
      if constexpr (std::is_invocable_v<Functor, ic<n0>>)
        f(ic<n0>{});
      else
        f();

      dpl::sfor<n0 + 1, n1>(f);
    }
  }

  template <int n, typename Functor>
  constexpr void sfor(const Functor& f) {
    dpl::sfor<0, n>(f);
  }

  namespace internal
  {
    template <int n0, int n1, typename Functor>
    constexpr void psfor(std::future<void>* queue, const Functor& f) {
      if constexpr (n0 < n1) {
        if constexpr (std::is_invocable_v<Functor, ic<n0>>)
          *queue = std::async(std::launch::async, f, ic<n0>{});
        else
          *queue = std::async(std::launch::async, f);        
          
        internal::psfor<n0 + 1, n1>(++queue, f);
      }
    }
  }
  
  template <int n, typename Functor>
  constexpr void psfor(const Functor& f) {
    std::array<std::future<void>, n> queue;

    internal::psfor<0, n>(queue.data(), f);

    sfor<n>([&](auto i) {
      queue[i].wait();
    });
    
    // for (std::future<void>& task : queue)
    //   task.wait();
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

  namespace extrapolant
  {
    static inline constexpr struct linear_t {} linear;
    static inline constexpr struct flat_t {} flat;
  }

  namespace lerp_base
  {
    struct linear {};
    struct log10 {};
  };
  
  
  
  template <typename Base = lerp_base::linear, typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(const T& a, const T& b, double t) {
    if constexpr (std::is_same_v<Extrapolant, extrapolant::linear_t> && std::is_same_v<Base, lerp_base::linear>)
      return a + (b - a)*t;
    else if constexpr (std::is_same_v<Base, lerp_base::log10>)
      return std::pow(10., dpl::lerp<lerp_base::linear, Extrapolant>(std::log10(a), std::log10(b), t));
    else if constexpr (std::is_same_v<Extrapolant, extrapolant::flat_t>)
      return
        t < 0 ? a :
        t > 1 ? b :
        dpl::lerp<Base, extrapolant::linear_t>(a, b, t);
    else
      static_assert(always_false<T>, "wrong extrapolant or base");
  }

  template <typename Base = lerp_base::linear, typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(T* ptr, double t) {
    return dpl::lerp<Base, Extrapolant>(ptr[0], ptr[1], t);
  }

  template <typename Extrapolant = extrapolant::linear_t, typename T>
  constexpr auto lerp(bool log10, const T& a, const T& b, double t) {
    return log10
      ? dpl::lerp<lerp_base::log10, Extrapolant>(a, b, t)
      : dpl::lerp<lerp_base::linear, Extrapolant>(a, b, t);
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
  



  // template <typename Func>
  // void apply(const Func& f) {
  //   f(value_);
  //   Base::apply(f);
  // }

  // template <typename Action>
  // void apply(const Action&) {}

  // template <typename Action, typename Item, typename... Rest>
  // void apply(const Action& a, Item&& item, Rest&&... args) {
  //   a(std::forward<Item>(item));
  //   dpl::apply(a, std::forward<Rest>(args)...);
  // }
  
  // inline void AddItem(QBoxLayout* ) { }
  //
  // template <typename... Rest>
  // void AddItem(QBoxLayout* l, QWidget& widget, Rest&&... args) {
  //   l->addWidget(&widget);
  //   QLayoutBuilder::AddItem(l, std::forward<Rest>(args)...);
  // }
}


template <std::integral T, typename Tag>
struct std::incrementable_traits<dpl::strong_integer<T, Tag>> : std::incrementable_traits<T> {};

template <std::integral T, typename Tag>
struct std::hash<dpl::strong_integer<T, Tag>>
{
  std::size_t operator()(const dpl::strong_integer<T, Tag>& k) const {
    return std::hash<T>()(*k);
  }
};

// namespace std
// {
//   // template <integral T, typename Tag>
//   // struct incrementable_traits<dpl::strong_integer<T, Tag>> : incrementable_traits<T> {};
//
//   // template <integral T, typename Tag>
//   // struct hash<dpl::strong_integer<T, Tag>>
//   // {
//   //   size_t operator()(const dpl::strong_integer<T, Tag>& k) const {
//   //     return hash<T>()(*k);
//   //   }
//   // };
// }
