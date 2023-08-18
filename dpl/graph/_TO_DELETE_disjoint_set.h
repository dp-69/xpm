#pragma once
#include <vector>

namespace HW { namespace dynamic_connectivity
{
  using namespace std;



  template <class _Category, class _Ty, class _Diff = ptrdiff_t, class _Pointer = _Ty*, class _Reference = _Ty&>
  struct TODO_MY_iterator { // base type for iterator classes
    using iterator_category = _Category;
    using value_type        = _Ty;
    using difference_type   = _Diff;
    using pointer           = _Pointer;
    using reference         = _Reference;
  };
  
  template <class Integral>
  class recursive_range_iterator : public TODO_MY_iterator<forward_iterator_tag, Integral>
  {
    Integral* _nextMap;
    Integral _current;

    constexpr static Integral END_VALUE = numeric_limits<Integral>::max();

  public:
    recursive_range_iterator(Integral* const nextMap, Integral const current)
      : _nextMap(nextMap),
        _current(current) {}

    static recursive_range_iterator end() {
      return recursive_range_iterator<Integral>(nullptr, END_VALUE);
    }
    

    recursive_range_iterator& operator++() {
      auto prev = _current;
      _current = _nextMap[_current];
      if (_current == prev)
        _current = END_VALUE;

      return *this;
    }

    recursive_range_iterator& operator++(int) {
      recursive_range_iterator temp(*this);
      operator++();
      return temp;
    }

    bool operator==(const recursive_range_iterator& rhs) {
      return _current == rhs._current;
    }

    bool operator!=(const recursive_range_iterator& rhs) {
      return _current != rhs._current;
    }

    Integral operator*() {
      return _current;
    }

    Integral* operator->() {
      return &_current;
    }
  };


  template <class Integral>
  class recursive_range
  {
    Integral* _nextMap;
    Integral _first;

    typedef recursive_range_iterator<Integral> iterator;

  public:
    recursive_range(Integral* const nextMap, Integral const first) : _nextMap(nextMap), _first(first) {}

    iterator begin() { return iterator(_nextMap, _first); }

    iterator end() { return iterator::end(); }
  };

  class disjoint_set
  {           
    size_t _count;
    vector<size_t> _rank;            
    vector<size_t> _tail;
    vector<size_t> _next;
    vector<size_t> _parent; 
    
    // i, j - representatives. i != j
    void link_rep_set(size_t rep_i, size_t rep_j) {            
      if (_rank[rep_i] > _rank[rep_j]) {
        _parent[rep_j] = rep_i;
        
        _next[_tail[rep_i]] = rep_j;
        _tail[rep_i] = _tail[rep_j];
      }
      else {
        _parent[rep_i] = rep_j;

        _next[_tail[rep_j]] = rep_i;
        _tail[rep_j] = _tail[rep_i];

        if (_rank[rep_i] == _rank[rep_j]) 
          _rank[rep_j]++;
      }
    }

  public:                
    explicit disjoint_set(size_t count) {      
      _count = count;     

      _parent.resize(count);
      _next.resize(count);
      _tail.resize(count);
      _rank.resize(count);

      for (size_t i = 0; i < count; i++) {
        _parent[i] = i;
        _next[i] = i;
        _tail[i] = i;
        _rank[i] = 0;
      } 
    }    

    recursive_range<size_t> range_of_set(size_t i) {
      return recursive_range<size_t>(_next.data(), find_rep(i));
    }     

    bool is_same_set(size_t i, size_t j) {
      return find_rep(i) == find_rep(j);
    }

    size_t find_rep(size_t elem) {
      auto old = elem;

      auto ancestor = _parent[elem];
      while (ancestor != elem) {
        elem = ancestor;
        ancestor = _parent[elem];
      }

      elem = _parent[old];
      
      while (ancestor != elem) {
        _parent[old] = ancestor;
        old = elem;
        elem = _parent[old];
      }

      return ancestor;
    }

    size_t number_of_components() const {
      size_t n = 0;

      for (size_t i = 0; i < _count; ++i)
        if (_parent[i] == i)
          ++n;

      return n;
    }

    // i, j - any element
    void union_set(size_t i, size_t j) {                  
      i = find_rep(i);
      j = find_rep(j);

      if (i != j)
        link_rep_set(i, j);
    }    
  };
}}
