#pragma once

#include "avltree_algorithms_ext.hpp"

namespace boost
{
  namespace intrusive
  {
    template <typename NodeTraits>
    class aug_avltree_algorithms_ext : public avltree_algorithms_ext<NodeTraits>
    {
      using base = avltree_algorithms_ext<NodeTraits>;
      using bstree_algo = bstree_algorithms<NodeTraits>;
      using nt = NodeTraits;

    public:
      using node_traits = nt;

      using node = typename nt::node;
      using node_ptr = typename nt::node_ptr;
      using const_node_ptr = typename nt::const_node_ptr;
      using balance = typename nt::balance;
      using subsize = typename nt::subsize;

      static void init(node_ptr node) {
        nt::set_size(node, 0);
        base::init(node);
      }

      static void size_rotate_right(node_ptr x, node_ptr x_oldleft) {
        nt::set_size(x_oldleft, nt::get_size(x));

        const_node_ptr x_right = nt::get_right(x);
        const_node_ptr x_oldleft_right = nt::get_right(x_oldleft);

        nt::set_size(x,
          (x_right ? nt::get_size(x_right) : 0) +
          (x_oldleft_right ? nt::get_size(x_oldleft_right) : 0) +
          1);
      }

      static void size_rotate_left(node_ptr x, node_ptr x_oldright) {
        nt::set_size(x_oldright, nt::get_size(x));

        const_node_ptr x_left = nt::get_left(x);
        const_node_ptr x_oldright_left = nt::get_left(x_oldright);

        nt::set_size(x,
          (x_left ? nt::get_size(x_left) : 0) +
          (x_oldright_left ? nt::get_size(x_oldright_left) : 0) +
          1);
      }

      static void size_increment_after_insert(const_node_ptr header, node_ptr n) {
        nt::set_size(n, 1);
        for (n = nt::get_parent(n); n != header; n = nt::get_parent(n))
          nt::set_size(n, nt::get_size(n) + 1);
      }

      static node_ptr avl_rotate_left_right(const node_ptr a, const node_ptr a_oldleft, const node_ptr& hdr) {
        size_rotate_left(a_oldleft, nt::get_right(a_oldleft));
        size_rotate_right(a, nt::get_right(a_oldleft));
        return base::avl_rotate_left_right(a, a_oldleft, hdr);
      }

      static node_ptr avl_rotate_right_left(const node_ptr a, const node_ptr a_oldright, node_ptr hdr) {
        size_rotate_right(a_oldright, nt::get_left(a_oldright));
        size_rotate_left(a, nt::get_left(a_oldright));
        return base::avl_rotate_right_left(a, a_oldright, hdr);
      }

      static void avl_rotate_left(const node_ptr& x, const node_ptr& x_oldright, const node_ptr& hdr) {
        size_rotate_left(x, x_oldright);
        base::avl_rotate_left(x, x_oldright, hdr);
      }

      static void avl_rotate_right(node_ptr x, node_ptr x_oldleft, node_ptr hdr) {
        size_rotate_right(x, x_oldleft);
        base::avl_rotate_right(x, x_oldleft, hdr);
      }

      static node_ptr erase(node_ptr header, node_ptr z) {       
        typename bstree_algo::data_for_rebalance info;

        bstree_algo::erase(header, z, info);
        if (info.y != z) {
          nt::set_balance(info.y, nt::get_balance(z));
          nt::set_size(info.y, nt::get_size(z));
        }
        
        for (node_ptr n = info.x_parent; n != header; n = nt::get_parent(n))
          nt::set_size(n, nt::get_size(n) - 1);

        rebalance_after_erasure(header, info.x, info.x_parent);
        return z;
      }

      static void push_back(node_ptr header, node_ptr new_node) {
        bstree_algo::push_back(header, new_node);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
      }

      static void push_front(node_ptr header, node_ptr new_node) {
        bstree_algo::push_front(header, new_node);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
      }

      template <typename NodePtrCompare>
      static node_ptr insert_equal_lower_bound(node_ptr header, node_ptr new_node, NodePtrCompare comp) {
        bstree_algo::insert_equal_lower_bound(header, new_node, comp);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      template <typename NodePtrCompare>
      static node_ptr insert_equal_upper_bound(node_ptr header, node_ptr new_node, NodePtrCompare comp) {
        bstree_algo::insert_equal_upper_bound(header, new_node, comp);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      static node_ptr insert_before(node_ptr header, node_ptr pos, node_ptr new_node) {
        bstree_algo::insert_before(header, pos, new_node);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }      





      static void join_trees(node_ptr hdr_left, node_ptr x, node_ptr hdr_right) {
        auto root_a = nt::get_parent(hdr_left);
        auto root_b = nt::get_parent(hdr_right);

        if (!root_b) {
          push_back(hdr_left, x);
          return;
        }

        if (!root_a) {
          bstree_algorithms<nt>::swap_tree(hdr_left, hdr_right);            
          push_front(hdr_left, x);
          return;
        }

        auto height_a = base::node_height(root_a);
        auto height_b = base::node_height(root_b);            
          
        if (abs(height_b - height_a) <= 1) {                        
          nt::set_parent(hdr_left, x);
          nt::set_left(hdr_left, nt::get_left(hdr_left));
          nt::set_right(hdr_left, nt::get_right(hdr_right));

          nt::set_parent(x, hdr_left);

          nt::set_left(x, root_a);
          nt::set_right(x, root_b);

          nt::set_parent(root_a, x);
          nt::set_parent(root_b, x);

          nt::set_balance(x, 
            height_a < height_b ? nt::positive() :
            height_a > height_b ? nt::negative() :
            nt::zero());

          nt::set_size(x, nt::get_size(root_a) + nt::get_size(root_b) + 1);

          base::init_header(hdr_right);

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

            nt::set_size(v, nt::get_size(v) + nt::get_size(root_b) + 1);

            v = nt::get_right(v);
          }

          auto v_parent = nt::get_parent(v);

          nt::set_left(x, v);
          nt::set_parent(v, x);

          nt::set_right(x, root_b);        
          nt::set_parent(root_b, x);
    
          nt::set_right(v_parent, x);
          nt::set_parent(x, v_parent);    

          nt::set_right(hdr_left, nt::get_right(hdr_right));

          nt::set_balance(x, height_b == h ? nt::zero() : nt::negative());

          nt::set_size(x, nt::get_size(v) + nt::get_size(root_b) + 1);
    
          base::init_header(hdr_right);
        }
        else {            
          node_ptr v = root_b;
          auto h = height_b;
          while (h > height_a + 1) {
            if (nt::get_balance(v) == nt::positive())
              h -= 2;
            else
              --h;

            nt::set_size(v, nt::get_size(v) + nt::get_size(root_a) + 1);

            v = nt::get_left(v);
          }

          node_ptr v_parent = nt::get_parent(v);

          nt::set_right(x, v);
          nt::set_parent(v, x);

          nt::set_left(x, root_a);        
          nt::set_parent(root_a, x);
    
          nt::set_left(v_parent, x);
          nt::set_parent(x, v_parent);    

          nt::set_left(hdr_right, nt::get_left(hdr_left));

          nt::set_balance(x, height_a == h ? nt::zero() : nt::positive());

          nt::set_size(x, nt::get_size(v) + nt::get_size(root_a) + 1);

          base::init_header(hdr_left);
          bstree_algo::swap_tree(hdr_left, hdr_right);
        }

        rebalance_after_insertion_no_balance_assignment(hdr_left, x);
      }


      static void split_tree(node_ptr hdr_left, node_ptr k, node_ptr hdr_right) {
        if (nt::get_right(hdr_left) == k) {
          erase(hdr_left, k);
          return;
        }

        if (nt::get_left(hdr_left) == k) {
          erase(hdr_left, k);
          bstree_algorithms<nt>::swap_tree(hdr_left, hdr_right);
          return;
        }

        
        auto k_height = avltree_algorithms_ext<NodeTraits>::node_height(k);

        auto left_tail_node = nt::get_left(k);
        auto left_tail_height = k_height - (nt::get_balance(k) == nt::positive() ? 2 : 1);

        auto right_tail_node = nt::get_right(k);
        auto right_tail_height = k_height - (nt::get_balance(k) == nt::negative() ? 2 : 1);

        node_ptr left_rightmost = nullptr;
        node_ptr right_leftmost = nullptr;

        if (left_tail_node) {
          left_rightmost = left_tail_node;
          while (nt::get_right(left_rightmost))
            left_rightmost = nt::get_right(left_rightmost);
        }

        if (right_tail_node) {
          right_leftmost = right_tail_node;
          while (nt::get_left(right_leftmost))
            right_leftmost = nt::get_left(right_leftmost);
        }


        node_ptr node = k;
        node_ptr parent = nt::get_parent(k);
        node_ptr grand_parent = nt::get_parent(parent);

        auto height = k_height;

        while (parent != hdr_left) {
          if (nt::get_left(parent) == node) {
            if (!right_leftmost)
              right_leftmost = parent;

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
              
              nt::set_size(x,
                (l_node ? nt::get_size(l_node) : 0) +
                (r_node ? nt::get_size(r_node) : 0) +
                1);

              right_tail_node = x;
            }
            else {
              if (l_height == -1) {
                nt::set_left(x, nullptr);
                nt::set_right(x, nullptr);
                nt::set_balance(x, nt::zero());
                
                nt::set_size(x, 1);
                nt::set_size(r_node, nt::get_size(r_node) + 1);

                node_ptr l = r_node;

                while (nt::get_left(l)) {
                  l = nt::get_left(l);
                  nt::set_size(l, nt::get_size(l) + 1);
                }

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

                  nt::set_size(v, nt::get_size(v) + nt::get_size(l_node) + 1);

                  v = nt::get_left(v);
                }

                auto v_parent = nt::get_parent(v);

                nt::set_right(x, v);
                nt::set_parent(v, x);

                nt::set_left(x, l_node);
                if (l_node)
                  nt::set_parent(l_node, x);

                nt::set_left(v_parent, x);
                nt::set_parent(x, v_parent);

                nt::set_balance(x, l_height == h ? nt::zero() : nt::positive());

                nt::set_size(x, nt::get_size(v) + nt::get_size(l_node) + 1);
              }

              nt::set_parent(hdr_right, r_node);
              nt::set_parent(r_node, hdr_right);

              if (rebalance_after_insertion_no_balance_assignment(hdr_right, x))
                right_tail_height = r_height;
              else
                right_tail_height = r_height + 1;

              right_tail_node = nt::get_parent(hdr_right);
            }
          }
          else {
            if (!left_rightmost)
              left_rightmost = parent;

            auto l_height =
              nt::get_balance(parent) == nt::negative()
                ? ++++height - 1
                : nt::get_balance(parent) == nt::positive()
                    ? ++height - 2
                    : ++height - 1;

            const auto& r_height = left_tail_height;

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

              nt::set_size(x,
                (l_node ? nt::get_size(l_node) : 0) +
                (r_node ? nt::get_size(r_node) : 0) +
                1);

              left_tail_node = x;
            }
            else {
              if (r_height == -1) {
                nt::set_left(x, nullptr);
                nt::set_right(x, nullptr);
                nt::set_balance(x, nt::zero());

                nt::set_size(x, 1);
                
                auto r = l_node;

                nt::set_size(l_node, nt::get_size(l_node) + 1);

                while (nt::get_right(r)) {
                  r = nt::get_right(r);
                  nt::set_size(r, nt::get_size(r) + 1);
                }

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

                  nt::set_size(v, nt::get_size(v) + nt::get_size(r_node) + 1);

                  v = nt::get_right(v);
                }

                auto v_parent = nt::get_parent(v);

                nt::set_left(x, v);
                nt::set_parent(v, x);

                nt::set_right(x, r_node);
                if (r_node)
                  nt::set_parent(r_node, x);

                nt::set_right(v_parent, x);
                nt::set_parent(x, v_parent);

                nt::set_balance(x, r_height == h ? nt::zero() : nt::negative());

                nt::set_size(x, nt::get_size(v) + nt::get_size(r_node) + 1);
              }

              nt::set_parent(hdr_right, l_node);
              nt::set_parent(l_node, hdr_right);

              if (rebalance_after_insertion_no_balance_assignment(hdr_right, x))
                left_tail_height = l_height;
              else
                left_tail_height = l_height + 1;

              left_tail_node = nt::get_parent(hdr_right);
            }
          }

          node = parent;
          parent = grand_parent;
          grand_parent = nt::get_parent(grand_parent);
        }

        nt::set_parent(hdr_right, right_tail_node);
        nt::set_left(hdr_right, right_leftmost);
        nt::set_right(hdr_right, nt::get_right(hdr_left));

        nt::set_parent(right_tail_node, hdr_right);


        nt::set_parent(hdr_left, left_tail_node);
        //nt::set_left(headerA, nt::get_left(headerA));                    
        nt::set_right(hdr_left, left_rightmost);

        nt::set_parent(left_tail_node, hdr_left);

        init(k);
      }



      
                     
      

      
      














      

     

     

//      static void init_header(const node_ptr& header) {
//        bstree_algorithms<nt>::init_header(header);
//        
//        nt::init_header(header);
//
////        NodeTraits::set_balance(header, NodeTraits::zero());
////        NodeTraits::augment_init(header);
//      }



//      static node_ptr erase(const node_ptr& header, const node_ptr& z) {
//        typename bstree_algorithms<nt>::data_for_rebalance info;
//
//
//        bstree_algorithms<nt>::erase(header, z, info);
//        if (info.y != z) {
//          NodeTraits::set_balance(info.y, NodeTraits::get_balance(z));
//        }
//        //Rebalance avltree
//        rebalance_after_erasure(header, info.x, info.x_parent);
//        return z;
//      }

      
      

      //   private:

     



      





     

      







      

      static void rebalance_after_erasure(node_ptr header, node_ptr x, node_ptr x_parent) {
        for (node_ptr root = NodeTraits::get_parent(header)
             ; x != root; root = NodeTraits::get_parent(header) , x_parent = NodeTraits::get_parent(x)) {
          const balance x_parent_balance = NodeTraits::get_balance(x_parent);
          //Don't cache x_is_leftchild or similar because x can be null and
          //equal to both x_parent_left and x_parent_right
          const node_ptr x_parent_left(NodeTraits::get_left(x_parent));
          const node_ptr x_parent_right(NodeTraits::get_right(x_parent));

          if (x_parent_balance == NodeTraits::zero()) {
            NodeTraits::set_balance(x_parent, x == x_parent_right ? NodeTraits::negative() : NodeTraits::positive());
            break; // the height didn't change, let's stop here
          }
          else if (x_parent_balance == NodeTraits::negative()) {
            if (x == x_parent_left) { ////x is left child or x and sibling are null
              NodeTraits::set_balance(x_parent, NodeTraits::zero()); // balanced
              x = x_parent;
            }
            else {
              // x is right child (x_parent_left is the left child)
              BOOST_INTRUSIVE_INVARIANT_ASSERT(x_parent_left);
              if (NodeTraits::get_balance(x_parent_left) == NodeTraits::positive()) {
                // x_parent_left MUST have a right child
                BOOST_INTRUSIVE_INVARIANT_ASSERT(NodeTraits::get_right(x_parent_left));
                x = avl_rotate_left_right(x_parent, x_parent_left, header);
              }
              else {
                avl_rotate_right(x_parent, x_parent_left, header);
                x = x_parent_left;
              }

              // if changed from negative to NodeTraits::positive(), no need to check above
              if (NodeTraits::get_balance(x) == NodeTraits::positive()) {
                break;
              }
            }
          }
          else if (x_parent_balance == NodeTraits::positive()) {
            if (x == x_parent_right) { //x is right child or x and sibling are null
              NodeTraits::set_balance(x_parent, NodeTraits::zero()); // balanced
              x = x_parent;
            }
            else {
              // x is left child (x_parent_right is the right child)
              BOOST_INTRUSIVE_INVARIANT_ASSERT(x_parent_right);
              if (NodeTraits::get_balance(x_parent_right) == NodeTraits::negative()) {
                // x_parent_right MUST have then a left child
                BOOST_INTRUSIVE_INVARIANT_ASSERT(NodeTraits::get_left(x_parent_right));
                x = avl_rotate_right_left(x_parent, x_parent_right, header);
              }
              else {
                avl_rotate_left(x_parent, x_parent_right, header);
                x = x_parent_right;
              }
              // if changed from NodeTraits::positive() to negative, no need to check above
              if (NodeTraits::get_balance(x) == NodeTraits::negative()) {
                break;
              }
            }
          }
          else {
            BOOST_INTRUSIVE_INVARIANT_ASSERT(false); // never reached
          }
        }
      }

      static bool rebalance_after_insertion(node_ptr header, node_ptr x) {
        NodeTraits::set_balance(x, NodeTraits::zero());
        return rebalance_after_insertion_no_balance_assignment(header, x);
      }

      static bool rebalance_after_insertion_no_balance_assignment(node_ptr header, node_ptr x) {
        // Rebalance.
        for (node_ptr root = NodeTraits::get_parent(header); x != root; root = NodeTraits::get_parent(header)) {
          node_ptr const x_parent(NodeTraits::get_parent(x));
          node_ptr const x_parent_left(NodeTraits::get_left(x_parent));
          const balance x_parent_balance = NodeTraits::get_balance(x_parent);
          const bool x_is_leftchild(x == x_parent_left);
          if (x_parent_balance == NodeTraits::zero()) {
            // if x is left, parent will have parent->bal_factor = negative
            // else, parent->bal_factor = NodeTraits::positive()
            NodeTraits::set_balance(x_parent, x_is_leftchild ? NodeTraits::negative() : NodeTraits::positive());
            x = x_parent;
          }
          else if (x_parent_balance == NodeTraits::positive()) {
            // if x is a left child, parent->bal_factor = zero
            if (x_is_leftchild)
              NodeTraits::set_balance(x_parent, NodeTraits::zero());
            else { // x is a right child, needs rebalancing
              if (NodeTraits::get_balance(x) == NodeTraits::negative())
                avl_rotate_right_left(x_parent, x, header);
              else
                avl_rotate_left(x_parent, x, header);
            }
            return true;
          }
          else if (x_parent_balance == NodeTraits::negative()) {
            // if x is a left child, needs rebalancing
            if (x_is_leftchild) {
              if (NodeTraits::get_balance(x) == NodeTraits::positive())
                avl_rotate_left_right(x_parent, x, header);
              else
                avl_rotate_right(x_parent, x, header);
            }
            else
              NodeTraits::set_balance(x_parent, NodeTraits::zero());
            return true;
          }
          else {
            BOOST_INTRUSIVE_INVARIANT_ASSERT(false); // never reached
          }
        }
        return false;
      }

      

      
    };  
  }
}







// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  DO NOT DELETE ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
//
//  static void refresh_augmentation(node_ptr n) {
//         if (n) {          
//           refresh_augmentation(nt::get_left(n));        
//           refresh_augmentation(nt::get_right(n));  
//
//           nt::augment_propagate(n, nt::get_left(n), nt::get_right(n));
//         }
//       }
//
//       static void refresh_size_recursive(const node_ptr& header) {
// //        refresh_size_recursive_rec(nt::get_parent(header));
//         refresh_augmentation(nt::get_parent(header));
//       }
//
//       static void refresh_size(const node_ptr& header) {               
//         refresh_size_iterative(header);
//       }
//
//       static void refresh_size_iterative(const node_ptr& header) {               
//         auto n = nt::get_parent(header);
//
//         if (!n)
//           return;
//
//         while (true) {
//           if (nt::get_left(n))
//             n = nt::get_left(n);
//           else if (nt::get_right(n))
//             n = nt::get_right(n);
//           else {              
// //            nt::set_size(n, 1);
//             nt::augment_identity(n);
//
//
//             auto p = nt::get_parent(n);                                        
//
//             while(true) {
//               if (p == header)
//                 return;              
//
//               if (!nt::get_right(p)) // no right child means n is the only left child
// //                nt::set_size(p, nt::get_size(n) + 1);
//                 nt::augment_propagate(p, n);
//               else if (nt::get_right(p) == n)
// //                nt::set_size(p, nt::get_size(nt::get_left(p)) + nt::get_size(n) + 1);
//                 nt::augment_propagate(p, nt::get_left(p), n);
//               else
//                 break;
//
//               n = p;
//               p = nt::get_parent(p);
//             }                         
//
//             n = nt::get_right(p);              
//           }         
//         }
//       }














// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  DO NOT DELETE ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
//  static subsize os_rank(const node_ptr& header, node_ptr x) {               
//         auto root = nt::get_parent(header);
//         auto idx = nt::get_size(nt::get_left(x));
//
//         while (x != root) {                    
//           auto parent = nt::get_parent(x);
//
//           if (x == nt::get_right(parent))
//             idx += nt::get_size(nt::get_left(parent)) + 1;
//
//           x = parent;
//         }
//
//         return idx;
//       }
//
//       static subsize os_rank_by_path(node_ptr* iter, node_ptr* path) {               
//         subsize idx = 0;
//                 
//         auto parent = *iter;
//         while (--iter != path) {
//           auto x = *iter;
//
//           if (x == nt::get_right(parent))
//             idx += nt::get_size(nt::get_left(parent)) + 1;
//
//
// //          auto exactRank = os_rank(header, *iter);
// //          auto currRank = idx + nt::get_size(nt::get_left(*iter));
// //
// //          if (exactRank != currRank) {
// //            auto err = 1;
// //          }
//
//           parent = x;
//         }
//
//         idx += nt::get_size(nt::get_left(*(path + 1))); 
//
//         return idx;
//       }
//
//       static node_ptr os_select(const node_ptr& header, subsize i) {
//         auto x = nt::get_parent(header);
//         
//         while (true) {
//           auto r = nt::get_size(nt::get_left(x));
//           if (i == r)
//             return x;
//           
//           if (i < r)
//             x = nt::get_left(x);
//           else {
//             x = nt::get_right(x);
//             i -= r + 1;
//           }
//         }
//       }