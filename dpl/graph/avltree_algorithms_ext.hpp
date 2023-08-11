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
      using nt = NodeTraits;

    public:
      using node_traits = NodeTraits;

      using node = typename nt::node;
      using node_ptr = typename nt::node_ptr;
      using const_node_ptr = typename nt::const_node_ptr;
      using balance = typename nt::balance;
      
      static const_node_ptr get_header(const_node_ptr node) {
        while (!is_header(node))
          node = nt::get_parent(node);
        return node;
      }

      static void init(node_ptr node) {
        base::init(node);
        nt::set_balance(node, nt::zero());
      }

      static void init_header(node_ptr header) {
        bstree_algorithms<nt>::init_header(header);        
        nt::init_header(header);
      }

      static bool is_header(const_node_ptr p) {
        return nt::is_header(p);
      }

      static auto node_height(const_node_ptr x) {
        int height = 0;

        while (nt::get_left(x) || nt::get_right(x)) {
          x = nt::get_balance(x) == nt::negative()
                ? nt::get_left(x)
                : nt::get_right(x);

          ++height;
        }

        return height;
      }

      static const_node_ptr* get_path(const_node_ptr n, const_node_ptr* path) {   
        while (!is_header(n)) {
          *++path = n; // TODO
          n = nt::get_parent(n);
        }       

        return path;
      }      
      
      /**
       * \brief l and r must not be equal
       */
      static bool less_than(const_node_ptr l, const_node_ptr r) {
        using default_path = default_path_buffer<nt>;

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
          return nt::get_right(*(r_subpath + 1)) == *r_subpath;

        return nt::get_left(*(l_subpath + 1)) == *l_subpath;        
      }    

       template<class UnaryLessThanNodeComparator>
      static node_ptr upper_bound(const_node_ptr header, UnaryLessThanNodeComparator comp) {
        node_ptr x = nt::get_parent(header);
        node_ptr y = header;

        while (x)
          if (comp(x)) {     // Key < x  
            y = x;
            x = nt::get_left(x);
          }
          else
            x = nt::get_right(x);

        return y;         
      }

      template<class UnaryMoreThanNodeComparator>
      static node_ptr lower_bound(const_node_ptr header, UnaryMoreThanNodeComparator comp) {
        const_node_ptr x = nt::get_parent(header);
        const_node_ptr y = header;

        while (x)
          if (comp(x))   // x < Key 
            x = nt::get_right(x);
          else {
            y = x;
            x = nt::get_left(x);
          }

        return y;         
      }

      static size_t calculate_subtree_size(const_node_ptr p) {
        return base::subtree_size(p);
      }

      static void join_trees(node_ptr header_a, node_ptr k, node_ptr header_b) {
        node_ptr root_a = nt::get_parent(header_a);
        node_ptr root_b = nt::get_parent(header_b);

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
          nt::set_parent(header_a, k);
          nt::set_left(header_a, nt::get_left(header_a));
          nt::set_right(header_a, nt::get_right(header_b));

          nt::set_parent(k, header_a);

          nt::set_left(k, root_a);
          nt::set_right(k, root_b);

          nt::set_parent(root_a, k);
          nt::set_parent(root_b, k);

          nt::set_balance(k,
            height_a < height_b ? nt::positive() :
            height_a > height_b ? nt::negative() :
            nt::zero()
          );

          init_header(header_b);
          return;            
        }
          

        if (height_a > height_b) {
          node_ptr v = root_a;
          auto h = height_a;
          while (h > height_b + 1) {
            if (nt::get_balance(v) == nt::negative())
              h -= 2;
            else
              --h;
            v = nt::get_right(v);
          }

          node_ptr v_parent = nt::get_parent(v);

          nt::set_left(k, v);
          nt::set_parent(v, k);

          nt::set_right(k, root_b);        
          nt::set_parent(root_b, k);

          nt::set_right(v_parent, k);
          nt::set_parent(k, v_parent);

          nt::set_balance(k, height_b == h ? nt::zero() : nt::negative());

          nt::set_right(header_a, nt::get_right(header_b));
    
          init_header(header_b);
        }
        else {            
          node_ptr v = root_b;
          auto h = height_b;
          while (h > height_a + 1) {
            if (nt::get_balance(v) == nt::positive())
              h -= 2;
            else
              --h;

            v = nt::get_left(v);
          }

          node_ptr v_parent = nt::get_parent(v);

          nt::set_right(k, v);
          nt::set_parent(v, k);

          nt::set_left(k, root_a);        
          nt::set_parent(root_a, k);
    
          nt::set_left(v_parent, k);
          nt::set_parent(k, v_parent);

          nt::set_balance(k, height_a == h ? nt::zero() : nt::positive());

          nt::set_left(header_b, nt::get_left(header_a));
    
          init_header(header_a);
          base::swap_tree(header_a, header_b);
        }

        base::rebalance_after_insertion_no_balance_assignment(header_a, k);
      }


      static void split_tree(node_ptr header_a, node_ptr k, node_ptr header_b) {
        if (nt::get_right(header_a) == k) {
          base::erase(header_a, k);
          return;
        }

        if (nt::get_left(header_a) == k) {
          base::erase(header_a, k);
          base::swap_tree(header_a, header_b);
          return;
        }

        auto k_height = node_height(k);

        node_ptr left_tail_node = nt::get_left(k);
        auto left_tail_height = k_height - (nt::get_balance(k) == nt::positive() ? 2 : 1);

        node_ptr right_tail_node = nt::get_right(k);
        auto right_tail_height = k_height - (nt::get_balance(k) == nt::negative() ? 2 : 1);

        node_ptr left_right = nullptr;
        node_ptr right_left = nullptr;

        if (left_tail_node) {
          left_right = left_tail_node;
          while (nt::get_right(left_right))
            left_right = nt::get_right(left_right);
        }

        if (right_tail_node) {
          right_left = right_tail_node;
          while (nt::get_left(right_left))
            right_left = nt::get_left(right_left);
        }


        node_ptr node = k;
        node_ptr parent = nt::get_parent(k);
        node_ptr grand_parent = nt::get_parent(parent);

        auto height = k_height;

        while (parent != header_a) {
          if (nt::get_left(parent) == node) {
            if (!right_left)
              right_left = parent;

            auto l_height = right_tail_height;
            auto r_height =
              nt::get_balance(parent) == nt::positive()
                ? ++++height - 1
                : nt::get_balance(parent) == nt::negative()
                    ? ++height - 2
                    : ++height - 1;

            node_ptr x = parent;
            node_ptr l_node = right_tail_node;
            node_ptr r_node = nt::get_right(x);

            auto height_diff = r_height - l_height;

            if (height_diff <= 1) {
              nt::set_left(x, l_node);
              if (l_node)
                nt::set_parent(l_node, x);

              if (height_diff == 1) {
                nt::set_balance(x, nt::positive());
                right_tail_height = r_height + 1;
              }
              else if (height_diff == 0) {
                nt::set_balance(x, nt::zero());
                right_tail_height = r_height + 1;
              }
              else /*if (heightDiff == -1)*/ {
                nt::set_balance(x, nt::negative());
                right_tail_height = r_height + 2;
              }

              right_tail_node = x;
            }
            else {
              if (l_height == -1) {
                nt::set_left(x, nullptr);
                nt::set_right(x, nullptr);
                nt::set_balance(x, nt::zero());
                
                node_ptr l = r_node;

                while (nt::get_left(l))
                  l = nt::get_left(l);

                nt::set_left(l, x);
                nt::set_parent(x, l);
              }
              else {
                node_ptr v = r_node;
                auto h = r_height;
                while (h > l_height + 1) {
                  if (nt::get_balance(v) == nt::positive())
                    h -= 2;
                  else
                    --h;                                    

                  v = nt::get_left(v);
                }

                node_ptr v_parent = nt::get_parent(v);

                nt::set_right(x, v);
                nt::set_parent(v, x);

                nt::set_left(x, l_node);
                if (l_node)
                  nt::set_parent(l_node, x);

                nt::set_left(v_parent, x);
                nt::set_parent(x, v_parent);

                nt::set_balance(x, l_height == h ? nt::zero() : nt::positive());
              }

              nt::set_parent(header_b, r_node);
              nt::set_parent(r_node, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                right_tail_height = r_height;
              else
                right_tail_height = r_height + 1;

              right_tail_node = nt::get_parent(header_b);
            }
          }
          else {
            if (!left_right)
              left_right = parent;

            auto l_height =
              nt::get_balance(parent) == nt::negative()
                ? ++++height - 1
                : nt::get_balance(parent) == nt::positive()
                    ? ++height - 2
                    : ++height - 1;

            auto r_height = left_tail_height;

            node_ptr x = parent;

            node_ptr l_node = nt::get_left(x);
            node_ptr r_node = left_tail_node;

            auto height_diff = l_height - r_height;

            if (height_diff <= 1) {
              nt::set_right(x, r_node);
              if (r_node)
                nt::set_parent(r_node, x);

              if (height_diff == 1) {
                nt::set_balance(x, nt::negative());
                left_tail_height = l_height + 1;
              }
              else if (height_diff == 0) {
                nt::set_balance(x, nt::zero());
                left_tail_height = l_height + 1;
              }
              else /*if (heightDiff == -1)*/ {
                nt::set_balance(x, nt::positive());
                left_tail_height = l_height + 2;
              }

              left_tail_node = x;
            }
            else {
              if (r_height == -1) {
                nt::set_left(x, nullptr);
                nt::set_right(x, nullptr);
                nt::set_balance(x, nt::zero());
                
                node_ptr r = l_node;

                while (nt::get_right(r))
                  r = nt::get_right(r);

                nt::set_right(r, x);
                nt::set_parent(x, r);               
              }
              else {
                node_ptr v = l_node;
                auto h = l_height;
                while (h > r_height + 1) {
                  if (nt::get_balance(v) == nt::negative())
                    h -= 2;
                  else
                    --h;

                  v = nt::get_right(v);
                }

                node_ptr v_parent = nt::get_parent(v);

                nt::set_left(x, v);
                nt::set_parent(v, x);

                nt::set_right(x, r_node);
                if (r_node)
                  nt::set_parent(r_node, x);

                nt::set_right(v_parent, x);
                nt::set_parent(x, v_parent);

                nt::set_balance(x, r_height == h ? nt::zero() : nt::negative());
              }

              nt::set_parent(header_b, l_node);
              nt::set_parent(l_node, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                left_tail_height = l_height;
              else
                left_tail_height = l_height + 1;

              left_tail_node = nt::get_parent(header_b);
            }
          }

          node = parent;
          parent = grand_parent;
          grand_parent = nt::get_parent(grand_parent);
        }

        nt::set_parent(header_b, right_tail_node);
        nt::set_left(header_b, right_left);
        nt::set_right(header_b, nt::get_right(header_a));

        nt::set_parent(right_tail_node, header_b);


        nt::set_parent(header_a, left_tail_node);
        /* nt::set_left(headerA, nt::get_left(headerA)); */                    
        nt::set_right(header_a, left_right);

        nt::set_parent(left_tail_node, header_a);

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
