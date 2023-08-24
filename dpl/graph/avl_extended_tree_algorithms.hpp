#pragma once

#include <boost/intrusive/set.hpp>

#include <boost/intrusive/avltree_algorithms.hpp>

namespace boost
{
  namespace intrusive
  {
    template<class NodeTraits>
    struct default_path_buffer
    {    
      typedef NodeTraits node_traits;
      typedef typename node_traits::node_ptr node_ptr;    

      static auto const ET_MAX_DEPTH = 256;

      static inline node_ptr path0[ET_MAX_DEPTH];
      static inline node_ptr path1[ET_MAX_DEPTH];
      static inline node_ptr path2[ET_MAX_DEPTH];
      static inline node_ptr path3[ET_MAX_DEPTH];
    };


    class avl_extended_tree_algorithms_not_implemented_error : public std::logic_error
    {
    public:
      avl_extended_tree_algorithms_not_implemented_error() : logic_error("not implemented") {}
    };






    template <class NodeTraits>
    class avl_extended_tree_algorithms : public avltree_algorithms<NodeTraits>
    {      
    public:
      using node_traits = NodeTraits;

      using node = typename node_traits::node;
      using node_ptr = node*;
      using const_node_ptr = const node*;
      using balance = typename NodeTraits::balance;

      using default_path = default_path_buffer<NodeTraits>;

      using base = avltree_algorithms<NodeTraits>;

      using bstree_algo = bstree_algorithms<NodeTraits>;
      using insert_commit_data = typename bstree_algorithms<NodeTraits>::insert_commit_data;     


      using this_algo = avl_extended_tree_algorithms<NodeTraits>;

      static node* get_header(node* node) {
        while (!is_header(node))
          node = node_traits::get_parent(node);
        return node;
      }

      static void init(const node_ptr& node) {
        bstree_algo::init(node);
        NodeTraits::set_balance(node, NodeTraits::zero());
      }

      static void init_header(const node_ptr& header) {
        bstree_algo::init_header(header);        
        node_traits::init_header(header);
//        NodeTraits::set_balance(header, NodeTraits::zero());        
      }

      static bool is_header(const node* p) {
        return node_traits::is_header(p);
      }


      static node_ptr* get_path(node_ptr node, node_ptr* path) {               
        while (!is_header(node)) {
          *++path = node;
          node = node_traits::get_parent(node);
        }       

        return path;
      }      
      
      static bool less_than(const node_ptr& n0, const node_ptr& n1) {
        auto iter0 = get_path(n0, default_path::path0);
        auto iter1 = get_path(n1, default_path::path1);
        return less_than(iter0, iter1, default_path::path0);
      }

      static bool less_than(node_ptr* iter0, node_ptr* iter1, node_ptr* path0) {                        
        // equal if *(path0 + 1) == *(path1 + 1)
        // should not be equal
        
        while (*iter0 == *iter1) {
          --iter0;
          --iter1;
        }

        if (iter0 == path0)
          return node_traits::get_right(*(iter1 + 1)) == *iter1;

        return node_traits::get_left(*(iter0 + 1)) == *iter0;        
      }    



      static void split_tree(const node_ptr& header_a, const node_ptr& split_node, const node_ptr& header_b) {
        if (node_traits::get_right(header_a) == split_node) {
          base::erase(header_a, split_node);
          return;
        }

        if (node_traits::get_left(header_a) == split_node) {
          base::erase(header_a, split_node);
          base::swap_tree(header_a, header_b);
          return;
        }

        auto splitNodeHeight = node_height(split_node);

        auto leftTailNode = node_traits::get_left(split_node);
        auto leftTailHeight = splitNodeHeight - (node_traits::get_balance(split_node) == node_traits::positive() ? 2 : 1);

        auto rightTailNode = node_traits::get_right(split_node);
        auto rightTailHeight = splitNodeHeight - (node_traits::get_balance(split_node) == node_traits::negative() ? 2 : 1);

        auto leftRight = node_ptr();
        auto rightLeft = node_ptr();

        if (leftTailNode) {
          leftRight = leftTailNode;
          while (node_traits::get_right(leftRight))
            leftRight = node_traits::get_right(leftRight);
        }

        if (rightTailNode) {
          rightLeft = rightTailNode;
          while (node_traits::get_left(rightLeft))
            rightLeft = node_traits::get_left(rightLeft);
        }


        auto node = split_node;
        auto parent = node_traits::get_parent(split_node);
        auto grandParent = node_traits::get_parent(parent);

        auto height = splitNodeHeight;

        while (parent != header_a) {
          if (node_traits::get_left(parent) == node) {
            if (!rightLeft)
              rightLeft = parent;

            const auto& lHeight = rightTailHeight;
            auto rHeight =
              node_traits::get_balance(parent) == node_traits::positive() ? ++++height - 1
                : node_traits::get_balance(parent) == node_traits::negative() ? ++height - 2 : ++height - 1;

            auto x = parent;

            const auto& lNode = rightTailNode;
            auto rNode = node_traits::get_right(x);

            auto heightDiff = rHeight - lHeight;

            if (heightDiff <= 1) {
              node_traits::set_left(x, lNode);
              if (lNode)
                node_traits::set_parent(lNode, x);

              if (heightDiff == 1) {
                node_traits::set_balance(x, node_traits::positive());
                rightTailHeight = rHeight + 1;
              }
              else if (heightDiff == 0) {
                node_traits::set_balance(x, node_traits::zero());
                rightTailHeight = rHeight + 1;
              }
              else /*if (heightDiff == -1)*/ {
                node_traits::set_balance(x, node_traits::negative());
                rightTailHeight = rHeight + 2;
              }

              
//              node_traits::set_size(x,
//                node_traits::get_size(lNode) + node_traits::get_size(rNode) + 1);

//              node_traits::augment_propagate(x, lNode, rNode);

              rightTailNode = x;
            }
            else {
              if (lHeight == -1) {
                node_traits::set_left(x, node_ptr());
                node_traits::set_right(x, node_ptr());
                node_traits::set_balance(x, node_traits::zero());
                
//                node_traits::set_size(x, 1);
//                node_traits::augment_identity(x);
                
                
                auto l = rNode;

//                node_traits::size_increment(rNode);
//                node_traits::augment_inserted_node(rNode, x);

                while (node_traits::get_left(l)) {
                  l = node_traits::get_left(l);
//                  node_traits::size_increment(l);
//                  node_traits::augment_inserted_node(l, x);
                }

                node_traits::set_left(l, x);
                node_traits::set_parent(x, l);
                
              }
              else {
//                auto lNodeSize = node_traits::get_size(lNode);

                auto v = rNode;
                auto hPrime = rHeight;
                while (hPrime > lHeight + 1) {
                  if (node_traits::get_balance(v) == node_traits::positive())
                    hPrime -= 2;
                  else
                    --hPrime;                                    

//                  node_traits::size_increment(v, lNodeSize + 1 /*x*/);

//                  node_traits::augment_inserted_node(v, x);
//                  node_traits::augment_inserted_subtree(v, lNode);                  

                  v = node_traits::get_left(v);
                }

                auto u = node_traits::get_parent(v);

                node_traits::set_right(x, v);
                node_traits::set_parent(v, x);

                node_traits::set_left(x, lNode);
                if (lNode)
                  node_traits::set_parent(lNode, x);

                node_traits::set_balance(x, lHeight == hPrime ? node_traits::zero() : node_traits::positive());

                node_traits::set_left(u, x);
                node_traits::set_parent(x, u);

//                node_traits::set_size(x, lNodeSize + node_traits::get_size(v) + 1);
//                node_traits::augment_propagate(x, lNode, v);
              }

              node_traits::set_parent(header_b, rNode);
              node_traits::set_parent(rNode, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                rightTailHeight = rHeight;
              else
                rightTailHeight = rHeight + 1;

              rightTailNode = node_traits::get_parent(header_b);
            }
          }
          else {
            if (!leftRight)
              leftRight = parent;

            auto lHeight =
              node_traits::get_balance(parent) == node_traits::negative() ? ++++height - 1
                : node_traits::get_balance(parent) == node_traits::positive() ? ++height - 2 : ++height - 1;

            const auto& rHeight = leftTailHeight;

            auto x = parent;

            auto lNode = node_traits::get_left(x);
            const auto& rNode = leftTailNode;

            auto heightDiff = lHeight - rHeight;

            if (heightDiff <= 1) {
              node_traits::set_right(x, rNode);
              if (rNode)
                node_traits::set_parent(rNode, x);

              if (heightDiff == 1) {
                node_traits::set_balance(x, node_traits::negative());
                leftTailHeight = lHeight + 1;
              }
              else if (heightDiff == 0) {
                node_traits::set_balance(x, node_traits::zero());
                leftTailHeight = lHeight + 1;
              }
              else /*if (heightDiff == -1)*/ {
                node_traits::set_balance(x, node_traits::positive());
                leftTailHeight = lHeight + 2;
              }

//              node_traits::set_size(x,
//                node_traits::get_size(lNode) + node_traits::get_size(rNode) + 1);

//              node_traits::augment_propagate(x, lNode, rNode);

              leftTailNode = x;
            }
            else {
              if (rHeight == -1) {
                node_traits::set_left(x, node_ptr());
                node_traits::set_right(x, node_ptr());
                node_traits::set_balance(x, node_traits::zero());

//                node_traits::set_size(x, 1);
//                node_traits::augment_identity(x);
                
                
                auto r = lNode;

//                node_traits::size_increment(lNode);
//                node_traits::augment_inserted_node(lNode, x);

                while (node_traits::get_right(r)) {
                  r = node_traits::get_right(r);
//                  node_traits::size_increment(r);
//                  node_traits::augment_inserted_node(r, x);
                }

                node_traits::set_right(r, x);
                node_traits::set_parent(x, r);               
              }
              else {
//                auto rNodeSize = node_traits::get_size(rNode);

                auto v = lNode;
                auto hPrime = lHeight;
                while (hPrime > rHeight + 1) {
                  if (node_traits::get_balance(v) == node_traits::negative())
                    hPrime -= 2;
                  else
                    --hPrime;

//                  node_traits::size_increment(v, rNodeSize + 1 /*x*/);

//                  node_traits::augment_inserted_node(v, x);
//                  node_traits::augment_inserted_subtree(v, rNode);

                  v = node_traits::get_right(v);
                }

                auto u = node_traits::get_parent(v);

                node_traits::set_left(x, v);
                node_traits::set_parent(v, x);

                node_traits::set_right(x, rNode);
                if (rNode)
                  node_traits::set_parent(rNode, x);

                node_traits::set_balance(x, rHeight == hPrime ? node_traits::zero() : node_traits::negative());

                node_traits::set_right(u, x);
                node_traits::set_parent(x, u);

//                node_traits::set_size(x, node_traits::get_size(v) + rNodeSize + 1);
//                node_traits::augment_propagate(x, v, rNode);
              }

              node_traits::set_parent(header_b, lNode);
              node_traits::set_parent(lNode, header_b);

              if (base::rebalance_after_insertion_no_balance_assignment(header_b, x))
                leftTailHeight = lHeight;
              else
                leftTailHeight = lHeight + 1;

              leftTailNode = node_traits::get_parent(header_b);
            }
          }

          node = parent;
          parent = grandParent;
          grandParent = node_traits::get_parent(grandParent);
        }

        node_traits::set_parent(header_b, rightTailNode);
        node_traits::set_left(header_b, rightLeft);
        node_traits::set_right(header_b, node_traits::get_right(header_a));

        node_traits::set_parent(rightTailNode, header_b);


        node_traits::set_parent(header_a, leftTailNode);
        //node_traits::set_left(headerA, node_traits::get_left(headerA));                    
        node_traits::set_right(header_a, leftRight);

        node_traits::set_parent(leftTailNode, header_a);

        init(split_node);
      }




      static void join_trees(const node_ptr& headerA, const node_ptr& x, const node_ptr& headerB) {
        auto rootA = node_traits::get_parent(headerA);
        auto rootB = node_traits::get_parent(headerB);

        if (!rootB) {
          base::push_back(headerA, x);
          return;
        }

        if (!rootA) {
          base::swap_tree(headerA, headerB);            
          base::push_front(headerA, x);
          return;
        }

        auto heightA = node_height(rootA);
        auto heightB = node_height(rootB);            
          
        if (abs(heightB - heightA) <= 1) {                        
          node_traits::set_parent(headerA, x);
          node_traits::set_left(headerA, node_traits::get_left(headerA));
          node_traits::set_right(headerA, node_traits::get_right(headerB));

          node_traits::set_parent(x, headerA);

          node_traits::set_left(x, rootA);
          node_traits::set_right(x, rootB);

          node_traits::set_parent(rootA, x);
          node_traits::set_parent(rootB, x);

          node_traits::set_balance(x,
            heightA < heightB ? node_traits::positive() :
            heightA > heightB ? node_traits::negative() :
            node_traits::zero()
          );

          // if (heightA < heightB)
          //   node_traits::set_balance(x, node_traits::positive());
          // else if (heightA > heightB)
          //   node_traits::set_balance(x, node_traits::negative());
          // else
          //   node_traits::set_balance(x, node_traits::zero());

//          node_traits::set_size(x,
//            node_traits::get_size(rootA) + node_traits::get_size(rootB) + 1);
//          node_traits::augment_propagate(x, rootA, rootB);

          init_header(headerB);
          return;            
        }
          



        if (heightA > heightB) {
//          auto sizeB = node_traits::get_size(rootB);
                        
          auto v = rootA;
          auto hPrime = heightA;
          while (hPrime > heightB + 1) {
            if (node_traits::get_balance(v) == node_traits::negative())
              hPrime -= 2;
            else
              --hPrime;

//            node_traits::size_increment(v, sizeB + 1 /*x*/);

//            node_traits::augment_inserted_node(v, x);
//            node_traits::augment_inserted_subtree(v, rootB);

            v = node_traits::get_right(v);
          }

          auto v_parent = node_traits::get_parent(v);

          node_traits::set_left(x, v);
          node_traits::set_parent(v, x);

          node_traits::set_right(x, rootB);        
          node_traits::set_parent(rootB, x);

          node_traits::set_right(v_parent, x);
          node_traits::set_parent(x, v_parent);

          node_traits::set_balance(x, heightB == hPrime ? node_traits::zero() : node_traits::negative());

          node_traits::set_right(headerA, node_traits::get_right(headerB));

//          node_traits::set_size(x, node_traits::get_size(v) + sizeB + 1);
//          node_traits::augment_propagate(x, v, rootB);
    
          init_header(headerB);
        }
        else {            
//          auto sizeA = node_traits::get_size(rootA);

          auto v = rootB;
          auto hPrime = heightB;
          while (hPrime > heightA + 1) {
            if (node_traits::get_balance(v) == node_traits::positive())
              hPrime -= 2;
            else
              --hPrime;

//            node_traits::size_increment(v, sizeA + 1 /*x*/);

//            node_traits::augment_inserted_node(v, x);
//            node_traits::augment_inserted_subtree(v, rootA);

            v = node_traits::get_left(v);
          }

          auto v_parent = node_traits::get_parent(v);

          node_traits::set_right(x, v);
          node_traits::set_parent(v, x);

          node_traits::set_left(x, rootA);        
          node_traits::set_parent(rootA, x);
    
          node_traits::set_left(v_parent, x);
          node_traits::set_parent(x, v_parent);

          node_traits::set_balance(x, heightA == hPrime ? node_traits::zero() : node_traits::positive());

          node_traits::set_left(headerB, node_traits::get_left(headerA));
    
//          node_traits::set_size(x, sizeA + node_traits::get_size(v) + 1);
//          node_traits::augment_propagate(x, rootA, v);

          
          init_header(headerA);
          base::swap_tree(headerA, headerB);
        }

        base::rebalance_after_insertion_no_balance_assignment(headerA, x);
      }

//      static bool join_trees(const node_ptr& headerA, const node_ptr& x, const node_ptr& headerB) {
//        auto rootA = node_traits::get_parent(headerA);
//        auto rootB = node_traits::get_parent(headerB);
//
//        if (!rootB)
//          return push_back(headerA, x);
//
//        if (!rootA) {
//          swap_tree(headerA, headerB);            
//          return push_front(headerA, x);
//        }
//
//        auto heightA = node_height(rootA);
//        auto heightB = node_height(rootB);            
//          
//        if (abs(heightB - heightA) <= 1) {                        
//          node_traits::set_parent(headerA, x);
//          node_traits::set_left(headerA, node_traits::get_left(headerA));
//          node_traits::set_right(headerA, node_traits::get_right(headerB));
//
//          node_traits::set_parent(x, headerA);
//
//          node_traits::set_left(x, rootA);
//          node_traits::set_right(x, rootB);
//
//          node_traits::set_parent(rootA, x);
//          node_traits::set_parent(rootB, x);
//   
//          if (heightA < heightB)
//            node_traits::set_balance(x, node_traits::positive());
//          else if (heightA > heightB)
//            node_traits::set_balance(x, node_traits::negative());
//          else
//            node_traits::set_balance(x, node_traits::zero());
//
////          node_traits::set_size(x,
////            node_traits::get_size(rootA) + node_traits::get_size(rootB) + 1);
////          node_traits::augment_propagate(x, rootA, rootB);
//
//          init_header(headerB);
//
//
//          return true;            
//        }
//          
//
//
//
//        if (heightA > heightB) {
////          auto sizeB = node_traits::get_size(rootB);
//                        
//          auto v = rootA;
//          auto hPrime = heightA;
//          while (hPrime > heightB + 1) {
//            if (node_traits::get_balance(v) == node_traits::negative())
//              hPrime -= 2;
//            else
//              --hPrime;
//
////            node_traits::size_increment(v, sizeB + 1 /*x*/);
//
////            node_traits::augment_inserted_node(v, x);
////            node_traits::augment_inserted_subtree(v, rootB);
//
//            v = node_traits::get_right(v);
//          }
//
//          auto u = node_traits::get_parent(v);
//
//          node_traits::set_left(x, v);
//          node_traits::set_parent(v, x);
//
//          node_traits::set_right(x, rootB);        
//          node_traits::set_parent(rootB, x);
//
//          node_traits::set_balance(x, heightB == hPrime ? node_traits::zero() : node_traits::negative());
//    
//          node_traits::set_right(u, x);
//          node_traits::set_parent(x, u);    
//
//          node_traits::set_right(headerA, node_traits::get_right(headerB));
//
////          node_traits::set_size(x, node_traits::get_size(v) + sizeB + 1);
////          node_traits::augment_propagate(x, v, rootB);
//    
//          init_header(headerB);
//        }
//        else {            
////          auto sizeA = node_traits::get_size(rootA);
//
//          auto v = rootB;
//          auto hPrime = heightB;
//          while (hPrime > heightA + 1) {
//            if (node_traits::get_balance(v) == node_traits::positive())
//              hPrime -= 2;
//            else
//              --hPrime;
//
////            node_traits::size_increment(v, sizeA + 1 /*x*/);
//
////            node_traits::augment_inserted_node(v, x);
////            node_traits::augment_inserted_subtree(v, rootA);
//
//            v = node_traits::get_left(v);
//          }
//
//          auto u = node_traits::get_parent(v);
//
//          node_traits::set_right(x, v);
//          node_traits::set_parent(v, x);
//
//          node_traits::set_left(x, rootA);        
//          node_traits::set_parent(rootA, x);
//
//          node_traits::set_balance(x, heightA == hPrime ? node_traits::zero() : node_traits::positive());
//    
//          node_traits::set_left(u, x);
//          node_traits::set_parent(x, u);    
//
//          node_traits::set_left(headerB, node_traits::get_left(headerA));
//    
////          node_traits::set_size(x, sizeA + node_traits::get_size(v) + 1);
////          node_traits::augment_propagate(x, rootA, v);
//
//          init_header(headerA);
//          swap_tree(headerA, headerB);
//        }
//
//        return rebalance_after_insertion_no_balance_assignment(headerA, x);
//      }
      
      template <class Integral = int>
      static Integral node_height(const node_ptr& node) {
        if (!node)
          return -1;

        Integral height = 0;

        auto x = node;          

        while (node_traits::get_left(x) || node_traits::get_right(x)) {
          x = node_traits::get_balance(x) == node_traits::negative() ? node_traits::get_left(x) : node_traits::get_right(x);
          ++height;
        }

        return height;
      }
               

     

      // DP
      template<class UnaryLessThanNodeComparator>
      static node_ptr upper_bound(const node_ptr& header, UnaryLessThanNodeComparator comp) {
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

      // DP
      //
      template<class UnaryMoreThanNodeComparator>
      static node_ptr lower_bound(const node_ptr& header, UnaryMoreThanNodeComparator comp) {
        node_ptr x = NodeTraits::get_parent(header);
        node_ptr y = header;

        while (x)
          if (comp(x))   // x < Key 
            x = NodeTraits::get_right(x);
          else {
            y = x;
            x = NodeTraits::get_left(x);
          }

        return y;         
      }









      


     

      

      static size_t calculate_subtree_size(const const_node_ptr& p) {
        return bstree_algo::subtree_size(p);
      }
    };
    


    const algo_types AvlExtendedTreeAlgorithms = algo_types(10000);

    template <class NodeTraits>
    struct get_algo<AvlExtendedTreeAlgorithms, NodeTraits>
    {
      typedef avl_extended_tree_algorithms<NodeTraits> type;
    };


    //    template <class ValueTraits, class NodePtrCompare, class ExtraChecker>
    //    struct get_node_checker<AvlAugmentedTreeAlgorithms, ValueTraits, NodePtrCompare, ExtraChecker>
    //    {
    //        typedef detail::avltree_node_checker<ValueTraits, NodePtrCompare, ExtraChecker> type;
    //    };
  }
}
