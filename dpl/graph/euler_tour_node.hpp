#pragma once

#include "avl_extended_augmented_tree_algorithms.hpp"

#include "cyclic_operations.hpp"
#include "dynamic_connectivity_graph.hpp"
#include "smart_pool.hpp"

namespace HW { namespace dynamic_connectivity
{   
  struct euler_tour_node : non_copyable_movable
  {     
//    string _TEXT_TO_REMOVE_;

  private:        
    typedef euler_tour_node* node_ptr;

//    friend struct smart_pool_traits<euler_tour_node>;
    friend struct euler_tour_node_traits;
    friend struct default_avl_lrpb_node_traits<euler_tour_node>;

    
    // has to be the first field for the smart_pool
    // size_t has to be 8 bytes = 64 bits, x64 system
    // [61 bits, element pointer] + [1 bit, entry type] + [2 bits for avl balance]        
    size_t ptr_type_balance_; 

    node_ptr parent_, left_, right_;
    


    static const auto type_bitmap_ = static_cast<size_t>(1u) << 2;
    static const auto balance_bitmap_ = type_bitmap_ - 1;
    static const auto ptr_bitmap_ = ~((static_cast<size_t>(1u) << 3) - 1);           
  };


  struct euler_tour_non_tree_edge_node;
      
  struct euler_tour_node_traits : default_avl_lrpb_node_traits<euler_tour_node>
  {    
  private:
    template<class T> 
    static T get_ptr(const node_ptr& n) {
      return reinterpret_cast<T>(n->ptr_type_balance_ & node::ptr_bitmap_);
    }

  public:
    typedef euler_tour_non_tree_edge_node* etnte_node_ptr;

    
    static etnte_node_ptr get_non_tree_edge_header(const node_ptr& n) {
      return get_ptr<etnte_node_ptr>(n);
    }

    static void set_non_tree_edge_header(const node_ptr& n, etnte_node_ptr etnteHeader) {
      n->ptr_type_balance_ = reinterpret_cast<size_t>(etnteHeader) | node::balance_bitmap_;
    }



    static bool is_loop_edge(const node_ptr& n) {
      return n->ptr_type_balance_ & node::type_bitmap_;
    }
    
    static directed_edge_ptr get_directed_edge(const node_ptr& n) {
      return get_ptr<directed_edge_ptr>(n);
    }

    static void set_directed_edge(const node_ptr& n, directed_edge_ptr directed_edge) {                  
      n->ptr_type_balance_ = reinterpret_cast<size_t>(directed_edge) | n->ptr_type_balance_ & node::balance_bitmap_;
    }

    static vertex_ptr get_vertex(const node_ptr& n) {
      return get_ptr<vertex_ptr>(n);
    } 
    
    static void set_vertex(const node_ptr& n, vertex_ptr vertex) {
      n->ptr_type_balance_ = reinterpret_cast<size_t>(vertex) | node::type_bitmap_ | n->ptr_type_balance_ & node::balance_bitmap_;      
    }                     
        
    static void init_header(const node_ptr& n) {
      n->ptr_type_balance_ = n->ptr_type_balance_ & ~node::balance_bitmap_ | node::balance_bitmap_;
    }

    // Has to be very efficient O(1) predicate.
    static bool is_header(const const_node_ptr& n) {
      return (n->ptr_type_balance_ & node::balance_bitmap_) == node::balance_bitmap_;
    }

    static balance get_balance(const const_node_ptr& n) {      
      return static_cast<balance>(n->ptr_type_balance_ & node::balance_bitmap_);
    }

    static void set_balance(const node_ptr& n, balance b) {
      n->ptr_type_balance_ = n->ptr_type_balance_ & ~node::balance_bitmap_ | static_cast<size_t>(b);
    }   
  };

      
  typedef bi::avl_extended_tree_algorithms<euler_tour_node_traits> euler_tour_algorithms;


  typedef cyclic_operations<euler_tour_algorithms> euler_tour_cyclic_operations;
  
 
//  class euler_tour_iterator : public iterator<forward_iterator_tag, euler_tour_node>
//  {    
//    typedef euler_tour_node_traits::node_ptr node_ptr;
//
//    node_ptr _node;      
//
//  public:
//    explicit euler_tour_iterator(const node_ptr& node)
//      : _node(node) {}
//
//    euler_tour_iterator& operator++() {
//      _node = euler_tour_cyclic_operations::next_node(_node);        
//      return *this;
//    }
//
//    euler_tour_iterator& operator++(int) {
//      auto temp(*this);
//      operator++();
//      return temp;
//    }
//
//    bool operator==(const euler_tour_iterator& rhs) const {
//      return _node == rhs._node;
//    }
//
//    bool operator!=(const euler_tour_iterator& rhs) const {
//      return _node != rhs._node;
//    }
//
//    node_ptr ptr() const {
//      return _node;
//    }
//
//    euler_tour_node& operator*() const {
//      return *_node;
//    }
//
//    node_ptr operator->() const {
//      return _node;
//    }
//  };




  typedef euler_tour_node et_node;
  typedef euler_tour_node_ptr et_node_ptr;
  typedef euler_tour_algorithms et_algo;
//  typedef bi::avl_extended_tree_algorithms<euler_tour_node_traits> et_algo_not_augmented;
  typedef euler_tour_node_traits et_nt;
  typedef et_nt::const_node_ptr et_const_node_ptr;
  typedef euler_tour_cyclic_operations et_cyclic_op;



   //  template<>
//  struct smart_pool_traits<et_node>
//  {
//    static et_node_ptr get_next(et_node_ptr x) {            
//      return static_cast<et_node_ptr>(x->graph_element_entry_type_avl_balance);
//    }
//
//    static void set_next(et_node_ptr x, et_node_ptr y) {
//      x->graph_element_entry_type_avl_balance = y;
//    }
//  };
}}
