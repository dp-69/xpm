#pragma once

#include <boost/intrusive/set.hpp>

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

      inline static node_ptr path0[ET_MAX_DEPTH];
      inline static node_ptr path1[ET_MAX_DEPTH];
      inline static node_ptr path2[ET_MAX_DEPTH];
      inline static node_ptr path3[ET_MAX_DEPTH];
    };




    class avl_extended_tree_algorithms_not_implemented_error : public std::logic_error
    {
    public:
      avl_extended_tree_algorithms_not_implemented_error() : logic_error("Function not yet implemented") { };
    };






    template <class NodeTraits>
    class avl_extended_tree_algorithms : public bstree_algorithms<NodeTraits>
    {      
    public:
      typedef typename NodeTraits::node node;
      typedef NodeTraits node_traits;
      typedef typename NodeTraits::node_ptr node_ptr;
      typedef typename NodeTraits::const_node_ptr const_node_ptr;
      typedef typename NodeTraits::balance balance;      

      typedef default_path_buffer<NodeTraits> default_path;

      typedef bstree_algorithms<NodeTraits> bstree_algo;
      typedef typename bstree_algorithms<NodeTraits>::insert_commit_data insert_commit_data;     


      using this_algo = avl_extended_tree_algorithms<NodeTraits>;

      static node_ptr get_header(node_ptr node) {
        while (!is_header(node))
          node = node_traits::get_parent(node);
        return node;
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



      static void split_tree(const node_ptr& headerA, const node_ptr& splitNode, const node_ptr& headerB) {
        if (node_traits::get_right(headerA) == splitNode) {
          erase(headerA, splitNode);
          return;
        }

        if (node_traits::get_left(headerA) == splitNode) {
          erase(headerA, splitNode);
          bstree_algorithms<NodeTraits>::swap_tree(headerA, headerB);
          return;
        }

        auto splitNodeHeight = node_height(splitNode);

        auto leftTailNode = node_traits::get_left(splitNode);
        auto leftTailHeight = splitNodeHeight - (node_traits::get_balance(splitNode) == node_traits::positive() ? 2 : 1);

        auto rightTailNode = node_traits::get_right(splitNode);
        auto rightTailHeight = splitNodeHeight - (node_traits::get_balance(splitNode) == node_traits::negative() ? 2 : 1);

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


        auto node = splitNode;
        auto parent = node_traits::get_parent(splitNode);
        auto grandParent = node_traits::get_parent(parent);

        auto height = splitNodeHeight;

        while (parent != headerA) {
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

              node_traits::set_parent(headerB, rNode);
              node_traits::set_parent(rNode, headerB);

              if (rebalance_after_insertion_no_balance_assignment(headerB, x))
                rightTailHeight = rHeight + 1;
              else
                rightTailHeight = rHeight;

              rightTailNode = node_traits::get_parent(headerB);
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

              node_traits::set_parent(headerB, lNode);
              node_traits::set_parent(lNode, headerB);

              if (rebalance_after_insertion_no_balance_assignment(headerB, x))
                leftTailHeight = lHeight + 1;
              else
                leftTailHeight = lHeight;

              leftTailNode = node_traits::get_parent(headerB);
            }
          }

          node = parent;
          parent = grandParent;
          grandParent = node_traits::get_parent(grandParent);
        }

        node_traits::set_parent(headerB, rightTailNode);
        node_traits::set_left(headerB, rightLeft);
        node_traits::set_right(headerB, node_traits::get_right(headerA));

        node_traits::set_parent(rightTailNode, headerB);


        node_traits::set_parent(headerA, leftTailNode);
        //node_traits::set_left(headerA, node_traits::get_left(headerA));                    
        node_traits::set_right(headerA, leftRight);

        node_traits::set_parent(leftTailNode, headerA);

        init(splitNode);
      }




      static bool join_trees(const node_ptr& headerA, const node_ptr& x, const node_ptr& headerB) {
        auto rootA = node_traits::get_parent(headerA);
        auto rootB = node_traits::get_parent(headerB);

        if (!rootB)
          return push_back(headerA, x);

        if (!rootA) {
          bstree_algorithms<NodeTraits>::swap_tree(headerA, headerB);            
          return push_front(headerA, x);
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
   
          if (heightA < heightB)
            node_traits::set_balance(x, node_traits::positive());
          else if (heightA > heightB)
            node_traits::set_balance(x, node_traits::negative());
          else
            node_traits::set_balance(x, node_traits::zero());

//          node_traits::set_size(x,
//            node_traits::get_size(rootA) + node_traits::get_size(rootB) + 1);
//          node_traits::augment_propagate(x, rootA, rootB);

          init_header(headerB);


          return true;            
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

          auto u = node_traits::get_parent(v);

          node_traits::set_left(x, v);
          node_traits::set_parent(v, x);

          node_traits::set_right(x, rootB);        
          node_traits::set_parent(rootB, x);

          node_traits::set_balance(x, heightB == hPrime ? node_traits::zero() : node_traits::negative());
    
          node_traits::set_right(u, x);
          node_traits::set_parent(x, u);    

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

          auto u = node_traits::get_parent(v);

          node_traits::set_right(x, v);
          node_traits::set_parent(v, x);

          node_traits::set_left(x, rootA);        
          node_traits::set_parent(rootA, x);

          node_traits::set_balance(x, heightA == hPrime ? node_traits::zero() : node_traits::positive());
    
          node_traits::set_left(u, x);
          node_traits::set_parent(x, u);    

          node_traits::set_left(headerB, node_traits::get_left(headerA));
    
//          node_traits::set_size(x, sizeA + node_traits::get_size(v) + 1);
//          node_traits::augment_propagate(x, rootA, v);

          
          init_header(headerA);
          bstree_algorithms<NodeTraits>::swap_tree(headerA, headerB);
        }

        return rebalance_after_insertion_no_balance_assignment(headerA, x);
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
               

          














      static void swap_nodes(const node_ptr& node1, const node_ptr& node2) {
        throw avl_extended_tree_algorithms_not_implemented_error();


        if (node1 == node2)
          return;

        node_ptr header1(bstree_algo::get_header(node1)), header2(bstree_algo::get_header(node2));
        swap_nodes(node1, header1, node2, header2);
      }

      static void swap_nodes(const node_ptr& node1, const node_ptr& header1, const node_ptr& node2, const node_ptr& header2) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        if (node1 == node2)
          return;

        bstree_algo::swap_nodes(node1, header1, node2, header2);
        //Swap balance
        balance c = NodeTraits::get_balance(node1);
        NodeTraits::set_balance(node1, NodeTraits::get_balance(node2));
        NodeTraits::set_balance(node2, c);
      }

      static void replace_node(const node_ptr& node_to_be_replaced, const node_ptr& new_node) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        if (node_to_be_replaced == new_node)
          return;
        replace_node(node_to_be_replaced, bstree_algo::get_header(node_to_be_replaced), new_node);
      }

//      static void swap_tree(const node_ptr& header1, const node_ptr& header2) {
//        bstree_algo::swap_tree(header1, header2);
//      }

      static void replace_node(const node_ptr& node_to_be_replaced, const node_ptr& header, const node_ptr& new_node) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        bstree_algo::replace_node(node_to_be_replaced, header, new_node);
        NodeTraits::set_balance(new_node, NodeTraits::get_balance(node_to_be_replaced));
      }

      static void unlink(const node_ptr& node) {
        throw avl_extended_tree_algorithms_not_implemented_error();

//        node_ptr x = NodeTraits::get_parent(node);
//        if (x) {
//          while (!is_header(x))
//            x = NodeTraits::get_parent(x);
//          erase(x, node);
//        }
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

      static node_ptr erase(const node_ptr& header, const node_ptr& z) {
        typename bstree_algo::data_for_rebalance info;


        bstree_algo::erase(header, z, info);
        if (info.y != z) {
          NodeTraits::set_balance(info.y, NodeTraits::get_balance(z));
        }
        //Rebalance avltree
        rebalance_after_erasure(header, info.x, info.x_parent);
        return z;
      }

      template <class Cloner, class Disposer>
      static void clone
      (const const_node_ptr& source_header, const node_ptr& target_header, Cloner cloner, Disposer disposer) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        avltree_node_cloner<NodeTraits, Cloner> new_cloner(cloner);
        bstree_algo::clone(source_header, target_header, new_cloner, disposer);
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









      template <class NodePtrCompare>
      static node_ptr insert_equal_upper_bound
      (const node_ptr& header, const node_ptr& new_node, NodePtrCompare comp) {
        bstree_algo::insert_equal_upper_bound(header, new_node, comp);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      template <class NodePtrCompare>
      static node_ptr insert_equal_lower_bound
      (const node_ptr& header, const node_ptr& new_node, NodePtrCompare comp) {
//        throw avl_extended_tree_algorithms_not_implemented_error();

        bstree_algo::insert_equal_lower_bound(header, new_node, comp);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      template <class NodePtrCompare>
      static node_ptr insert_equal
      (const node_ptr& header, const node_ptr& hint, const node_ptr& new_node, NodePtrCompare comp) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        bstree_algo::insert_equal(header, hint, new_node, comp);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }

      static node_ptr insert_before
      (const node_ptr& header, const node_ptr& pos, const node_ptr& new_node) {
//        throw avl_extended_tree_algorithms_not_implemented_error();

        bstree_algo::insert_before(header, pos, new_node);
        rebalance_after_insertion(header, new_node);
        return new_node;
      }      

      static bool push_back(const node_ptr& header, const node_ptr& new_node) {
        bstree_algo::push_back(header, new_node);
        return rebalance_after_insertion(header, new_node);
      }      

      static bool push_front(const node_ptr& header, const node_ptr& new_node) {
        bstree_algo::push_front(header, new_node);       
        return rebalance_after_insertion(header, new_node);
      }

      static void insert_unique_commit
      (const node_ptr& header, const node_ptr& new_value, const insert_commit_data& commit_data) {
        throw avl_extended_tree_algorithms_not_implemented_error();

        bstree_algo::insert_unique_commit(header, new_value, commit_data);
        rebalance_after_insertion(header, new_value);
      }

      static bool is_header(const const_node_ptr& p) {
        return node_traits::is_header(p);
      }

      static size_t calculate_subtree_size(const const_node_ptr& p) {
        return bstree_algo::subtree_size(p);
      }

      static bool verify(const node_ptr& header) {
        std::size_t height;
        std::size_t count;
        return verify_recursion(NodeTraits::get_parent(header), count, height);
      }

      //   private:








      static bool verify_recursion(node_ptr n, std::size_t& count, std::size_t& height) {
        if (!n) {
          count = 0;
          height = 0;
          return true;
        }       

        std::size_t leftcount, rightcount;
        std::size_t leftheight, rightheight;
       

        if (!verify_recursion(NodeTraits::get_left(n), leftcount, leftheight) ||
            !verify_recursion(NodeTraits::get_right(n), rightcount, rightheight)) {
          return false;
        }
        count = 1u + leftcount + rightcount;
        height = 1u + (leftheight > rightheight ? leftheight : rightheight);

        //If equal height, balance must be zero
        if (rightheight == leftheight) {
          if (NodeTraits::get_balance(n) != NodeTraits::zero()) {
            BOOST_ASSERT(0);
            return false;
          }
        }
        //If right is taller than left, then the difference must be at least 1 and the balance positive
        else if (rightheight > leftheight) {
          if (rightheight - leftheight > 1) {
            BOOST_ASSERT(0);
            return false;
          }
          else if (NodeTraits::get_balance(n) != NodeTraits::positive()) {
            auto bal = NodeTraits::get_balance(n);
            
            BOOST_ASSERT(0);
            return false;
          }
        }
        //If left is taller than right, then the difference must be at least 1 and the balance negative
        else {
          if (leftheight - rightheight > 1) {
            BOOST_ASSERT(0);
            return false;
          }
          else if (NodeTraits::get_balance(n) != NodeTraits::negative()) {
            BOOST_ASSERT(0);
            return false;
          }
        }

        return true;
      }

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

      static bool rebalance_after_insertion(const node_ptr& header, node_ptr x) {
        NodeTraits::set_balance(x, NodeTraits::zero());
        return rebalance_after_insertion_no_balance_assignment(header, x);
      }

      static void left_right_balancing(const node_ptr& a, const node_ptr& b, const node_ptr& c) {
        // balancing...
        const balance c_balance = NodeTraits::get_balance(c);
        const balance zero_balance = NodeTraits::zero();
        const balance posi_balance = NodeTraits::positive();
        const balance nega_balance = NodeTraits::negative();
        NodeTraits::set_balance(c, zero_balance);
        if (c_balance == nega_balance) {
          NodeTraits::set_balance(a, posi_balance);
          NodeTraits::set_balance(b, zero_balance);
        }
        else if (c_balance == zero_balance) {
          NodeTraits::set_balance(a, zero_balance);
          NodeTraits::set_balance(b, zero_balance);
        }
        else if (c_balance == posi_balance) {
          NodeTraits::set_balance(a, zero_balance);
          NodeTraits::set_balance(b, nega_balance);
        }
        else {
          BOOST_INTRUSIVE_INVARIANT_ASSERT(false); // never reached
        }
      }

      static node_ptr avl_rotate_left_right(const node_ptr a, const node_ptr a_oldleft, const node_ptr& hdr) { // [note: 'a_oldleft' is 'b']
        //             |                               |         //
        //             a(-2)                           c         //
        //            / \                             / \        //
        //           /   \        ==>                /   \       //
        //      (pos)b    [g]                       b     a      //
        //          / \                            / \   / \     //
        //        [d]  c                         [d]  e f  [g]   //
        //            / \                                        //
        //           e   f                                       //
        const node_ptr c = NodeTraits::get_right(a_oldleft);

        bstree_algo::rotate_left_no_parent_fix(a_oldleft, c);
        //No need to link c with a [NodeTraits::set_parent(c, a) + NodeTraits::set_left(a, c)]
        //as c is not root and another rotation is coming
        
        bstree_algo::rotate_right(a, c, NodeTraits::get_parent(a), hdr);
        left_right_balancing(a, a_oldleft, c);
        return c;
      }

      static node_ptr avl_rotate_right_left(const node_ptr a, const node_ptr a_oldright, const node_ptr& hdr) { // [note: 'a_oldright' is 'b']
        //              |                               |           //
        //              a(pos)                          c           //
        //             / \                             / \          //
        //            /   \                           /   \         //
        //          [d]   b(neg)         ==>         a     b        //
        //               / \                        / \   / \       //
        //              c  [g]                    [d] e  f  [g]     //
        //             / \                                          //
        //            e   f                                         //
        const node_ptr c(NodeTraits::get_left(a_oldright));
        
        bstree_algo::rotate_right_no_parent_fix(a_oldright, c);
        //No need to link c with a [NodeTraits::set_parent(c, a) + NodeTraits::set_right(a, c)]
        //as c is not root and another rotation is coming.

        bstree_algo::rotate_left(a, c, NodeTraits::get_parent(a), hdr);
        left_right_balancing(a_oldright, a, c);
        return c;
      }

      static void avl_rotate_left(const node_ptr& x, const node_ptr& x_oldright, const node_ptr& hdr) {        
        bstree_algo::rotate_left(x, x_oldright, NodeTraits::get_parent(x), hdr);

        // reset the balancing factor
        if (NodeTraits::get_balance(x_oldright) == NodeTraits::positive()) {
          NodeTraits::set_balance(x, NodeTraits::zero());
          NodeTraits::set_balance(x_oldright, NodeTraits::zero());
        }
        else { // this doesn't happen during insertions
          NodeTraits::set_balance(x, NodeTraits::positive());
          NodeTraits::set_balance(x_oldright, NodeTraits::negative());
        }
      }


      static void avl_rotate_right(const node_ptr& x, const node_ptr& x_oldleft, const node_ptr& hdr) {        
        bstree_algo::rotate_right(x, x_oldleft, NodeTraits::get_parent(x), hdr);

        // reset the balancing factor
        if (NodeTraits::get_balance(x_oldleft) == NodeTraits::negative()) {
          NodeTraits::set_balance(x, NodeTraits::zero());
          NodeTraits::set_balance(x_oldleft, NodeTraits::zero());
        }
        else { // this doesn't happen during insertions
          NodeTraits::set_balance(x, NodeTraits::negative());
          NodeTraits::set_balance(x_oldleft, NodeTraits::positive());
        }
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
