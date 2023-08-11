#pragma once

#include "avltree_algorithms_ext.hpp"

namespace boost
{
  namespace intrusive
  {
    template <class NodeTraits>
    class avl_extended_augmented_tree_algorithms : public avltree_algorithms_ext<NodeTraits>
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

        nt::set_size(x,
          nt::get_size(nt::get_right(x)) +
          nt::get_size(nt::get_right(x_oldleft)) +
          1);
      }

      static void size_rotate_left(node_ptr x, node_ptr x_oldright) {
        nt::set_size(x_oldright, nt::get_size(x));

        nt::set_size(x,
          nt::get_size(nt::get_left(x)) +
          nt::get_size(nt::get_left(x_oldright)) +
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

      template <class NodePtrCompare>
      static node_ptr insert_equal_lower_bound(node_ptr header, node_ptr new_node, NodePtrCompare comp) {
        bstree_algo::insert_equal_lower_bound(header, new_node, comp);
        size_increment_after_insert(header, new_node);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      template <class NodePtrCompare>
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










      static void split_tree(const node_ptr& headerA, const node_ptr& splitNode, const node_ptr& headerB) {
        if (nt::get_right(headerA) == splitNode) {
          erase(headerA, splitNode);
          return;
        }

        if (nt::get_left(headerA) == splitNode) {
          erase(headerA, splitNode);
          bstree_algorithms<nt>::swap_tree(headerA, headerB);
          return;
        }

        
        auto splitNodeHeight = avltree_algorithms_ext<NodeTraits>::node_height(splitNode);

        auto leftTailNode = nt::get_left(splitNode);
        auto leftTailHeight = splitNodeHeight - (nt::get_balance(splitNode) == nt::positive() ? 2 : 1);

        auto rightTailNode = nt::get_right(splitNode);
        auto rightTailHeight = splitNodeHeight - (nt::get_balance(splitNode) == nt::negative() ? 2 : 1);

        auto leftRight = node_ptr();
        auto rightLeft = node_ptr();

        if (leftTailNode) {
          leftRight = leftTailNode;
          while (nt::get_right(leftRight))
            leftRight = nt::get_right(leftRight);
        }

        if (rightTailNode) {
          rightLeft = rightTailNode;
          while (nt::get_left(rightLeft))
            rightLeft = nt::get_left(rightLeft);
        }


        auto node = splitNode;
        auto parent = nt::get_parent(splitNode);
        auto grandParent = nt::get_parent(parent);

        auto height = splitNodeHeight;

        while (parent != headerA) {
          if (nt::get_left(parent) == node) {
            if (!rightLeft)
              rightLeft = parent;

            const auto& lHeight = rightTailHeight;
            auto rHeight =
              nt::get_balance(parent) == nt::positive() ? ++++height - 1
                : nt::get_balance(parent) == nt::negative() ? ++height - 2 : ++height - 1;

            auto x = parent;

            const auto& lNode = rightTailNode;
            auto rNode = nt::get_right(x);

            auto heightDiff = rHeight - lHeight;

            if (heightDiff <= 1) {
              nt::set_left(x, lNode);
              if (lNode)
                nt::set_parent(lNode, x);

              if (heightDiff == 1) {
                nt::set_balance(x, nt::positive());
                rightTailHeight = rHeight + 1;
              }
              else if (heightDiff == 0) {
                nt::set_balance(x, nt::zero());
                rightTailHeight = rHeight + 1;
              }
              else /*if (heightDiff == -1)*/ {
                nt::set_balance(x, nt::negative());
                rightTailHeight = rHeight + 2;
              }

              
//              nt::set_size(x,
//                nt::get_size(lNode) + nt::get_size(rNode) + 1);

              nt::augment_propagate(x, lNode, rNode);

              rightTailNode = x;
            }
            else {
              if (lHeight == -1) {
                nt::set_left(x, node_ptr());
                nt::set_right(x, node_ptr());
                nt::set_balance(x, nt::zero());
                
//                nt::set_size(x, 1);
                nt::augment_identity(x);
                
                
                auto l = rNode;

//                nt::size_increment(rNode);
                nt::augment_inserted_node(rNode, x);

                while (nt::get_left(l)) {
                  l = nt::get_left(l);
//                  nt::size_increment(l);
                  nt::augment_inserted_node(l, x);
                }

                nt::set_left(l, x);
                nt::set_parent(x, l);
                
              }
              else {
//                auto lNodeSize = nt::get_size(lNode);

                auto v = rNode;
                auto hPrime = rHeight;
                while (hPrime > lHeight + 1) {
                  if (nt::get_balance(v) == nt::positive())
                    hPrime -= 2;
                  else
                    --hPrime;                                    

//                  nt::size_increment(v, lNodeSize + 1 /*x*/);

                  nt::augment_inserted_node(v, x);
                  nt::augment_inserted_subtree(v, lNode);                  

                  v = nt::get_left(v);
                }

                auto u = nt::get_parent(v);

                nt::set_right(x, v);
                nt::set_parent(v, x);

                nt::set_left(x, lNode);
                if (lNode)
                  nt::set_parent(lNode, x);

                nt::set_balance(x, lHeight == hPrime ? nt::zero() : nt::positive());

                nt::set_left(u, x);
                nt::set_parent(x, u);

//                nt::set_size(x, lNodeSize + nt::get_size(v) + 1);
                nt::augment_propagate(x, lNode, v);
              }

              nt::set_parent(headerB, rNode);
              nt::set_parent(rNode, headerB);

              if (rebalance_after_insertion_no_balance_assignment(headerB, x))
                rightTailHeight = rHeight + 1;
              else
                rightTailHeight = rHeight;

              rightTailNode = nt::get_parent(headerB);
            }
          }
          else {
            if (!leftRight)
              leftRight = parent;

            auto lHeight =
              nt::get_balance(parent) == nt::negative() ? ++++height - 1
                : nt::get_balance(parent) == nt::positive() ? ++height - 2 : ++height - 1;

            const auto& rHeight = leftTailHeight;

            auto x = parent;

            auto lNode = nt::get_left(x);
            const auto& rNode = leftTailNode;

            auto heightDiff = lHeight - rHeight;

            if (heightDiff <= 1) {
              nt::set_right(x, rNode);
              if (rNode)
                nt::set_parent(rNode, x);

              if (heightDiff == 1) {
                nt::set_balance(x, nt::negative());
                leftTailHeight = lHeight + 1;
              }
              else if (heightDiff == 0) {
                nt::set_balance(x, nt::zero());
                leftTailHeight = lHeight + 1;
              }
              else /*if (heightDiff == -1)*/ {
                nt::set_balance(x, nt::positive());
                leftTailHeight = lHeight + 2;
              }

//              nt::set_size(x,
//                nt::get_size(lNode) + nt::get_size(rNode) + 1);

              nt::augment_propagate(x, lNode, rNode);

              leftTailNode = x;
            }
            else {
              if (rHeight == -1) {
                nt::set_left(x, node_ptr());
                nt::set_right(x, node_ptr());
                nt::set_balance(x, nt::zero());

//                nt::set_size(x, 1);
                nt::augment_identity(x);
                
                
                auto r = lNode;

//                nt::size_increment(lNode);
                nt::augment_inserted_node(lNode, x);

                while (nt::get_right(r)) {
                  r = nt::get_right(r);
//                  nt::size_increment(r);
                  nt::augment_inserted_node(r, x);
                }

                nt::set_right(r, x);
                nt::set_parent(x, r);               
              }
              else {
//                auto rNodeSize = nt::get_size(rNode);

                auto v = lNode;
                auto hPrime = lHeight;
                while (hPrime > rHeight + 1) {
                  if (nt::get_balance(v) == nt::negative())
                    hPrime -= 2;
                  else
                    --hPrime;

//                  nt::size_increment(v, rNodeSize + 1 /*x*/);

                  nt::augment_inserted_node(v, x);
                  nt::augment_inserted_subtree(v, rNode);

                  v = nt::get_right(v);
                }

                auto u = nt::get_parent(v);

                nt::set_left(x, v);
                nt::set_parent(v, x);

                nt::set_right(x, rNode);
                if (rNode)
                  nt::set_parent(rNode, x);

                nt::set_balance(x, rHeight == hPrime ? nt::zero() : nt::negative());

                nt::set_right(u, x);
                nt::set_parent(x, u);

//                nt::set_size(x, nt::get_size(v) + rNodeSize + 1);
                nt::augment_propagate(x, v, rNode);
              }

              nt::set_parent(headerB, lNode);
              nt::set_parent(lNode, headerB);

              if (rebalance_after_insertion_no_balance_assignment(headerB, x))
                leftTailHeight = lHeight + 1;
              else
                leftTailHeight = lHeight;

              leftTailNode = nt::get_parent(headerB);
            }
          }

          node = parent;
          parent = grandParent;
          grandParent = nt::get_parent(grandParent);
        }

        nt::set_parent(headerB, rightTailNode);
        nt::set_left(headerB, rightLeft);
        nt::set_right(headerB, nt::get_right(headerA));

        nt::set_parent(rightTailNode, headerB);


        nt::set_parent(headerA, leftTailNode);
        //nt::set_left(headerA, nt::get_left(headerA));                    
        nt::set_right(headerA, leftRight);

        nt::set_parent(leftTailNode, headerA);

        init(splitNode);
      }



      static void join_trees(const node_ptr& headerA, const node_ptr& x, const node_ptr& headerB) {
        auto rootA = nt::get_parent(headerA);
        auto rootB = nt::get_parent(headerB);

        if (!rootB) {
          push_back(headerA, x);
          return;
        }

        if (!rootA) {
          bstree_algorithms<nt>::swap_tree(headerA, headerB);            
          push_front(headerA, x);
          return;
        }

        auto heightA = avltree_algorithms_ext<NodeTraits>::node_height(rootA);
        auto heightB = avltree_algorithms_ext<NodeTraits>::node_height(rootB);            
          
        if (abs(heightB - heightA) <= 1) {                        
          nt::set_parent(headerA, x);
          nt::set_left(headerA, nt::get_left(headerA));
          nt::set_right(headerA, nt::get_right(headerB));

          nt::set_parent(x, headerA);

          nt::set_left(x, rootA);
          nt::set_right(x, rootB);

          nt::set_parent(rootA, x);
          nt::set_parent(rootB, x);
   
          if (heightA < heightB)
            nt::set_balance(x, nt::positive());
          else if (heightA > heightB)
            nt::set_balance(x, nt::negative());
          else
            nt::set_balance(x, nt::zero());

//          nt::set_size(x,
//            nt::get_size(rootA) + nt::get_size(rootB) + 1);
          nt::augment_propagate(x, rootA, rootB);

          avltree_algorithms_ext<NodeTraits>::init_header(headerB);

          return;            
        }
          



        if (heightA > heightB) {
//          auto sizeB = nt::get_size(rootB);
                        
          auto v = rootA;
          auto hPrime = heightA;
          while (hPrime > heightB + 1) {
            if (nt::get_balance(v) == nt::negative())
              hPrime -= 2;
            else
              --hPrime;

//            nt::size_increment(v, sizeB + 1 /*x*/);

            nt::augment_inserted_node(v, x);
            nt::augment_inserted_subtree(v, rootB);

            v = nt::get_right(v);
          }

          auto u = nt::get_parent(v);

          nt::set_left(x, v);
          nt::set_parent(v, x);

          nt::set_right(x, rootB);        
          nt::set_parent(rootB, x);

          nt::set_balance(x, heightB == hPrime ? nt::zero() : nt::negative());
    
          nt::set_right(u, x);
          nt::set_parent(x, u);    

          nt::set_right(headerA, nt::get_right(headerB));

//          nt::set_size(x, nt::get_size(v) + sizeB + 1);
          nt::augment_propagate(x, v, rootB);
    
          avltree_algorithms_ext<NodeTraits>::init_header(headerB);
        }
        else {            
//          auto sizeA = nt::get_size(rootA);

          auto v = rootB;
          auto hPrime = heightB;
          while (hPrime > heightA + 1) {
            if (nt::get_balance(v) == nt::positive())
              hPrime -= 2;
            else
              --hPrime;

//            nt::size_increment(v, sizeA + 1 /*x*/);

            nt::augment_inserted_node(v, x);
            nt::augment_inserted_subtree(v, rootA);

            v = nt::get_left(v);
          }

          auto u = nt::get_parent(v);

          nt::set_right(x, v);
          nt::set_parent(v, x);

          nt::set_left(x, rootA);        
          nt::set_parent(rootA, x);

          nt::set_balance(x, heightA == hPrime ? nt::zero() : nt::positive());
    
          nt::set_left(u, x);
          nt::set_parent(x, u);    

          nt::set_left(headerB, nt::get_left(headerA));
    
//          nt::set_size(x, sizeA + nt::get_size(v) + 1);
          nt::augment_propagate(x, rootA, v);

          avltree_algorithms_ext<NodeTraits>::init_header(headerA);
          bstree_algorithms<nt>::swap_tree(headerA, headerB);
        }

        rebalance_after_insertion_no_balance_assignment(headerA, x);
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

     



      





     

      







      

      static void rebalance_after_erasure(const node_ptr& header, node_ptr x, node_ptr x_parent) {
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

      static bool rebalance_after_insertion(const node_ptr& header, node_ptr x) {
        NodeTraits::set_balance(x, NodeTraits::zero());
        return rebalance_after_insertion_no_balance_assignment(header, x);
      }

      static bool rebalance_after_insertion_no_balance_assignment(const node_ptr& header, node_ptr x) {
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
            return false;
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
            return false;
          }
          else {
            BOOST_INTRUSIVE_INVARIANT_ASSERT(false); // never reached
          }
        }
        return true;
      }

      

      
    };  




    // const algo_types AvlExtendedAugmentedTreeAlgorithms = algo_types(10001);
    //
    // template <class NodeTraits>
    // struct get_algo<AvlExtendedAugmentedTreeAlgorithms, NodeTraits>
    // {
    //   typedef avl_extended_augmented_tree_algorithms<NodeTraits> type;
    // };




    //    template <class ValueTraits, class NodePtrCompare, class ExtraChecker>
    //    struct get_node_checker<AvlAugmentedTreeAlgorithms, ValueTraits, NodePtrCompare, ExtraChecker>
    //    {
    //        typedef detail::avltree_node_checker<ValueTraits, NodePtrCompare, ExtraChecker> type;
    //    };
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