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
// #include <list>
#include <memory>
#include <set>





namespace dpl
{
 

  

  
  
  template <typename ... Args>
  class callback
  {
    
    
  public:
    const int priority_;
    
  explicit callback(int p) : priority_(p) {}

  virtual ~callback() = default;

    virtual void call(Args&& ... args) = 0;


    bool operator< (const callback& r) const {
      return priority_ < r.priority_;
    }
    
    // bool operator<(const callback& rhs) const { return priority_ < rhs.priority_; }
    
    // friend bool operator<(callback& lhs, callback& rhs) { return lhs.priority_ < rhs.priority_; }
  // friend bool operator<=(const callback& lhs, const callback& rhs) { return !(rhs < lhs); }
  // friend bool operator>(const callback& lhs, const callback& rhs) { return rhs < lhs; }
  // friend bool operator>=(const callback& lhs, const callback& rhs) { return !(lhs < rhs); }
  // virtual void operator()() = 0;
  };


  
  

  template <typename Functor, typename ... Args>
  class callback_generic final : public callback<Args...>
  {
    Functor func_;
  public:
    explicit callback_generic(Functor func, int priority)
      : callback<Args...>(priority), func_(func) {}

    void call(Args&& ... args) override {
      if constexpr(std::is_invocable_v<Functor, Args...>)
        func_(std::forward<Args>(args)...);
      else
        func_();
    }
  };

  
  
  template <typename ... Args>
  class event
  {
    using _callback = callback<Args...>;
    using _unique_ptr_callback = std::unique_ptr<_callback>;
    
    struct callback_comparer
    {
      constexpr bool operator()(const _unique_ptr_callback& lhs, const _unique_ptr_callback& rhs) const {
        return *lhs < *rhs;
      }
    };

    
    std::multiset<std::unique_ptr<callback<Args...>>, callback_comparer> registered_;

  public:
    template <typename Callback>
    void add(std::unique_ptr<Callback> c) {
      registered_.insert(std::move(c));
    }
    
    template <typename Functor>
    void add(Functor func, int priority = 0) {
      this->add(std::make_unique<callback_generic<Functor, Args...>>(func, priority));
    }

    template <typename T>
    event& operator+=(T&& t){
      this->add(std::forward<T>(t));
      return *this;
    }
    
    void notify(Args&& ... args) {
      for (auto& callback : registered_)
        callback->call(std::forward<Args>(args)...);
    }

    void operator()(Args&& ... args) {
      this->notify(std::forward<Args>(args)...);
    }
  };
}


