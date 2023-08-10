#pragma once

#include <boost/intrusive/set.hpp>
#include <boost/intrusive/avltree_algorithms.hpp>

namespace boost
{
  namespace intrusive
  {
    template<class NodeTraits>
    struct default_path_buffer // TODO: maybe move to the function iself, check performance
    {
      using node_ptr = typename NodeTraits::node_ptr;
      using const_node_ptr = typename NodeTraits::const_node_ptr;

      static inline constexpr auto ET_MAX_DEPTH = 256;

      static inline const_node_ptr path0[ET_MAX_DEPTH];
      static inline const_node_ptr path1[ET_MAX_DEPTH];
      static inline const_node_ptr path2[ET_MAX_DEPTH];
      static inline const_node_ptr path3[ET_MAX_DEPTH];
    };


    template <typename NodeTraits>
    class avltree_algorithms_ext : public avltree_algorithms<NodeTraits>
    {
      using base = avltree_algorithms<NodeTraits>;

    public:
      using node_traits = NodeTraits;

      using node = typename node_traits::node;
      using node_ptr = typename node_traits::node_ptr;
      using const_node_ptr = typename node_traits::const_node_ptr;
      using balance = typename node_traits::balance;
      
      static const_node_ptr get_header(const_node_ptr node) {
        while (!is_header(node))
          node = node_traits::get_parent(node);
        return node;
      }

      static void init(node_ptr node) {
        base::init(node);
        node_traits::set_balance(node, node_traits::zero());
      }

      static void init_header(node_ptr header) {
        bstree_algorithms<node_traits>::init_header(header);        
        node_traits::init_header(header);
      }

      static bool is_header(const_node_ptr p) {
        return node_traits::is_header(p);
      }

      static auto node_height(const_node_ptr x) {
        int height = 0;

        while (node_traits::get_left(x) || node_traits::get_right(x)) {
          x = node_traits::get_balance(x) == node_traits::negative()
                ? node_traits::get_left(x)
                : node_traits::get_right(x);

          ++height;
        }

        return height;
      }

      static const_node_ptr* get_path(const_node_ptr n, const_node_ptr* path) {   
        while (!is_header(n)) {
          *++path = n; // TODO
          n = node_traits::get_parent(n);
        }       

        return path;
      }      
      
      /**
       * \brief l and r must not be equal
       */
      static bool less_than(const_node_ptr l, const_node_ptr r) {
        using default_path = default_path_buffer<node_traits>;

        auto l_subpath = get_path(l, default_path::path0);
        auto r_subpath = get_path(r, default_path::path1);
        return less_than(l_subpath, r_subpath, default_path::path0);
      }

      static bool less_than(const_node_ptr* l_subpath, const_node_ptr* r_subpath, const_node_ptr* l_path) {  // TODO                      
        // equal if *(path0 + 1) == *(path1 + 1) // TODO - try to understand
        // should not be equal
        
        while (*l_subpath == *r_subpath) {
          --l_subpath;
          --r_subpath;
        }

        if (l_subpath == l_path)
          return node_traits::get_right(*(r_subpath + 1)) == *r_subpath;

        return node_traits::get_left(*(l_subpath + 1)) == *l_subpath;        
      }    

       template<class UnaryLessThanNodeComparator>
      static node_ptr upper_bound(const_node_ptr header, UnaryLessThanNodeComparator comp) {
        node_ptr x = NodeTraits::get_parent(header);
        node_ptr y = header;

        while (x)
          if (comp(x)) {     // Key < x  
            y = x;
            x = NodeTraits::get_left(x);
          }
          else
            x = NodeTraits::get_right(x);

        return y;         
      }

      template<class UnaryMoreThanNodeComparator>
      static node_ptr lower_bound(const_node_ptr header, UnaryMoreThanNodeComparator comp) {
        const_node_ptr x = NodeTraits::get_parent(header);
        const_node_ptr y = header;

        while (x)
          if (comp(x))   // x < Key 
            x = NodeTraits::get_right(x);
          else {
            y = x;
            x = NodeTraits::get_left(x);
          }

        return y;         
      }

      static size_t calculate_subtree_size(const_node_ptr p) {
        return base::subtree_size(p);
      }

      static void join_trees(node_ptr header_a, node_ptr k, node_ptr header_b) {
        node_ptr root_a = node_traits::get_parent(header_a);
        node_ptr root_b = node_traits::get_parent(header_b);

        if (!root_b) {
          base::push_back(header_a, k);
          return;
        }

        if (!root_a) {
          base::swap_tree(header_a, header_b);            
          base::push_front(header_a, k);
          return;
        }

        auto height_a = node_height(root_a);
        auto height_b = node_height(root_b);            
          
        if (std::abs(height_b - height_a) <= 1) {                        
          node_traits::set_parent(header_a, k);
          node_traits::set_left(header_a, node_traits::get_left(header_a));
          node_traits::set_right(header_a, node_traits::get_right(header_b));

          node_traits::set_parent(k, header_a);

          node_traits::set_left(k, root_a);
          node_traits::set_right(k, root_b);

          node_traits::set_parent(root_a, k);
          node_traits::set_parent(root_b, k);

          node_traits::set_balance(k,
            height_a < height_b ? node_traits::positive() :
            height_a > height_b ? node_traits::negative() :
            node_traits::zero()
          );

          init_header(header_b);
          return;            
        }
          

        if (height_a > height_b) {
          node_ptr v = root_a;
          auto h = height_a;
          while (h > height_b + 1) {
            if (node_traits::get_balance(v) == node_traits::negative())
              h -= 2;
            else
              --h;
            v = node_traits::get_right(v);
          }

          node_ptr v_parent = node_traits::get_parent(v);

          node_traits::set_left(k, v);
          node_traits::set_parent(v, k);

          node_traits::set_right(k, root_b);        
          node_traits::set_parent(root_b, k);

          node_traits::set_right(v_parent, k);
          node_traits::set_parent(k, v_parent);

          node_traits::set_balance(k, height_b == h ? node_traits::zero() : node_traits::negative());

          node_traits::set_right(header_a, node_traits::get_right(header_b));
    
          init_header(header_b);
        }
        else {            
          node_ptr v = root_b;
          auto h = height_b;
          while (h > height_a + 1) {
            if (node_traits::get_balance(v) == node_traits::positive())
              h -= 2;
            else
              --h;

            v = node_traits::get_left(v);
          }

          node_ptr v_parent = node_traits::get_parent(v);

          node_traits::set_right(k, v);
          node_traits::set_parent(v, k);

          node_traits::set_left(k, root_a);        
          node_traits::set_parent(root_a, k);
    
          node_traits::set_left(v_parent, k);
          node_traits::set_parent(k, v_parent);

          node_traits::set_balance(k, height_a == h ? node_traits::zero() : node_traits::positive());

          node_traits::set_left(header_b, node_traits::get_left(header_a));
    
          init_header(header_a);
          base::swap_tree(header_a, header_b);
        }

        base::rebalance_after_insertion_no_balance_assignment(header_a, k);
      }


      static void split_tree(node_ptr header_a, node_ptr k, node_ptr header_b) {
        if (node_traits::get_right(header_a) == k) {
          base::erase(header_a, k);
          return;
        }

        if (node_traits::get_left(header_a) == k) {
          base::erase(header_a, k);
          base::swap_tree(header_a, header_b);
          return;
        }

        auto k_height = node_height(k);

        node_ptr left_tail_node = node_traits::get_left(k);
        auto left_tail_height = k_height - (node_traits::get_balance(k) == node_traits::positive() ? 2 : 1);

        node_ptr right_tail_node = node_traits::get_right(k);
        auto right_tail_height = k_height - (node_traits::get_balance(k) == node_traits::negative() ? 2 : 1);

        node_ptr left_right = nullptr;
        node_ptr right_left = nullptr;

        if (left_tail_node) {
          left_right = left_tail_node;
          while (node_traits::get_right(left_right))
            left_right = node_traits::get_right(left_right);
        }

        if (right_tail_node) {
          right_left = right_tail_node;
          while (node_traits::get_left(right_left))
            right_left = node_traits::get_left(right_left);
        }


        node_ptr node = k;
        node_ptr parent = node_traits::get_parent(k);
        node_ptr grand_parent = node_traits::get_parent(parent);

        auto height = k_height;

        while (parent != header_a) {
          if (node_traits::get_left(parent) == node) {
            if (!right_left)
              right_left = parent;

            auto l_height = right_tail_height;
            auto r_height =
              node_traits::get_balance(parent) == node_traits::positive()
                ? ++++height - 1
                : node_traits::get_balance(parent) == node_traits::negative()
                    ? ++height - 2
                    : ++height - 1;

            node_ptr x = parent;
            node_ptr l_node = right_tail_node;
            node_ptr r_node = node_traits::get_right(x);

            auto height_diff = r_height - l_height;

            if (height_diff <= 1) {
              node_traits::set_left(x, l_node);
              if (l_node)
                node_traits::set_parent(l_node, x);

              if (height_diff == 1) {
                node_traits::set_balance(x, node_traits::positive());
                right_tail_height = r_height + 1;
              }
              else if (height_diff == 0) {
                node_traits::set_balance(x, node_traits::zero());
                right_tail_height = r_height + 1;
              }
              else /*if (heightDiff == -1)*/ {
                node_traits::set_balance(x, node_traits::negative());
                right_tail_height = r_height + 2;
              }

              right_tail_node = x;
            }
            else {
              if (l_height == -1) {
                node_traits::set_left(x, nullptr);
                node_traits::set_right(x, nullptr);
                node_traits::set_balance(x, node_traits::zero());
                
                node_ptr l = r_node;

                while (node_traits::get_left(l))
                  l = node_traits::get_left(l);

                node_traits::set_left(l, x);
                node_traits::set_parent(x, l);
              }
              else {
                node_ptr v = r_node;
                auto h = r_height;
                while (h > l_height + 1) {
                  if (node_traits::get_balance(v) == node_traits::positive())
                    h -= 2;
                  else
                    --h;                                    

                  v = node_traits::get_left(v);
                }

                node_ptr v_parent = node_traits::get_parent(v);

                node_traits::set_right(x, v);
                node_traits::set_parent(v, x);

                node_traits::set_left(x, l_node);
                if (l_node)
                  node_traits::set_parent(l_node, x);

                node_traits::set_left(v_parent, x);
                node_traits::set_parent(x, v_parent);

                node_traits::set_balance(x, l_height == h ? node_traits::zero() : node_traits::positive());
              }

              node_traits::set_parent(header_b, r_node);
              node_traits::set_parent(r_node, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                right_tail_height = r_height;
              else
                right_tail_height = r_height + 1;

              right_tail_node = node_traits::get_parent(header_b);
            }
          }
          else {
            if (!left_right)
              left_right = parent;

            auto l_height =
              node_traits::get_balance(parent) == node_traits::negative()
                ? ++++height - 1
                : node_traits::get_balance(parent) == node_traits::positive()
                    ? ++height - 2
                    : ++height - 1;

            auto r_height = left_tail_height;

            node_ptr x = parent;

            node_ptr l_node = node_traits::get_left(x);
            node_ptr r_node = left_tail_node;

            auto height_diff = l_height - r_height;

            if (height_diff <= 1) {
              node_traits::set_right(x, r_node);
              if (r_node)
                node_traits::set_parent(r_node, x);

              if (height_diff == 1) {
                node_traits::set_balance(x, node_traits::negative());
                left_tail_height = l_height + 1;
              }
              else if (height_diff == 0) {
                node_traits::set_balance(x, node_traits::zero());
                left_tail_height = l_height + 1;
              }
              else /*if (heightDiff == -1)*/ {
                node_traits::set_balance(x, node_traits::positive());
                left_tail_height = l_height + 2;
              }

              left_tail_node = x;
            }
            else {
              if (r_height == -1) {
                node_traits::set_left(x, nullptr);
                node_traits::set_right(x, nullptr);
                node_traits::set_balance(x, node_traits::zero());
                
                node_ptr r = l_node;

                while (node_traits::get_right(r))
                  r = node_traits::get_right(r);

                node_traits::set_right(r, x);
                node_traits::set_parent(x, r);               
              }
              else {
                node_ptr v = l_node;
                auto h = l_height;
                while (h > r_height + 1) {
                  if (node_traits::get_balance(v) == node_traits::negative())
                    h -= 2;
                  else
                    --h;

                  v = node_traits::get_right(v);
                }

                node_ptr v_parent = node_traits::get_parent(v);

                node_traits::set_left(x, v);
                node_traits::set_parent(v, x);

                node_traits::set_right(x, r_node);
                if (r_node)
                  node_traits::set_parent(r_node, x);

                node_traits::set_right(v_parent, x);
                node_traits::set_parent(x, v_parent);

                node_traits::set_balance(x, r_height == h ? node_traits::zero() : node_traits::negative());
              }

              node_traits::set_parent(header_b, l_node);
              node_traits::set_parent(l_node, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                left_tail_height = l_height;
              else
                left_tail_height = l_height + 1;

              left_tail_node = node_traits::get_parent(header_b);
            }
          }

          node = parent;
          parent = grand_parent;
          grand_parent = node_traits::get_parent(grand_parent);
        }

        node_traits::set_parent(header_b, right_tail_node);
        node_traits::set_left(header_b, right_left);
        node_traits::set_right(header_b, node_traits::get_right(header_a));

        node_traits::set_parent(right_tail_node, header_b);


        node_traits::set_parent(header_a, left_tail_node);
        /* node_traits::set_left(headerA, node_traits::get_left(headerA)); */                    
        node_traits::set_right(header_a, left_right);

        node_traits::set_parent(left_tail_node, header_a);

        init(k);
      }




      


      
      
               

     

     
    };
    


    // const algo_types AvlExtendedTreeAlgorithms = algo_types(10000);
    //
    // template <class NodeTraits>
    // struct get_algo<AvlExtendedTreeAlgorithms, NodeTraits>
    // {
    //   typedef avl_extended_tree_algorithms<NodeTraits> type;
    // };
  }
}
