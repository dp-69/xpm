// #pragma once
//
// #include <boost/intrusive/list.hpp>
// #include <boost/intrusive/trivial_value_traits.hpp>
//
// namespace HW::dynamic_connectivity
// {
//   struct et_node;
//
//   using euler_tour_node_ptr = HW::dynamic_connectivity::et_node*;
//
//   struct dynamic_connectivity_graph;
// }
//
//
// namespace HW::dynamic_connectivity
// {
//   struct vertex;
//   using vertex_ptr = vertex*;
//
//   struct directed_edge;
//
//   namespace bi = boost::intrusive;
//
//
//
//   typedef bi::list<
//     directed_edge,
//     bi::value_traits<bi::trivial_value_traits<default_list_node_traits<directed_edge>, bi::normal_link>>,
//     bi::constant_time_size<false>> out_edge_list;
//
//   using out_edge_iterator = inorder_iter<default_list_node_traits<directed_edge>>;
//
//
//
//
//
//
//   struct vertex : non_copyable_movable
//   {
//     typedef vertex* node_ptr;
//
//     // intrusive list of graph vertices
//     node_ptr prev_;
//     node_ptr next_;          
//     
//     out_edge_list out_edges_;
//
//     
//
//   private:    
//     friend struct default_list_node_traits<vertex>;
//
//     
//
//     
//     friend void clear_out_edges(vertex_ptr, dynamic_connectivity_graph&);
//     friend void add_edge(vertex&, vertex&, directed_edge&, directed_edge&, dynamic_connectivity_graph&);
//     friend pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_ptr, const dynamic_connectivity_graph&);
//
//     friend out_edge_iterator out_edges_begin(const vertex_ptr);
//     friend out_edge_iterator out_edges_end(const vertex_ptr);
//
//
//
//   public:    
//     size_t row_idx_ = 0;  // relative index, not unique identifier    
//
//         
//     euler_tour_node_ptr et_entry_;
//     bool visited = false;
//   };
//
//
//
//   struct directed_edge : non_copyable_movable
//   {
//     typedef directed_edge* node_ptr;    
// //    typedef const directed_edge* const_node_ptr;
//
//     // entry for the list of out edges from v0 pointing to v1
//     node_ptr prev_;
//     node_ptr next_;
//
//   private:
//     friend struct default_list_node_traits<directed_edge>;
//
//     
//
//
//
// //    typedef pointer_plus_aligned_least_significant_bits<void*, bool, 1> compression;
//     typedef tagged_pointer_as_size_t<bool, 1> compression;
//     typedef euler_tour_node_ptr et_node_ptr;
//     typedef euler_tour_non_tree_edge_node_ptr etnte_node_ptr;
//
//     std::size_t entry_type_;
//
//   public:
//     vertex_ptr v1;  // pointing-in vertex  
//
//     node_ptr opposite;
//
//     
//     static void set_null_et_entry(const node_ptr& x) {
//       x->entry_type_ = 0;
//     }
//
//     static bool is_null_et_entry(const node_ptr& x) {
//       return !x->entry_type_;
//     }    
//
//
//
//     static bool is_tree_edge(const node_ptr& x) {
//       return compression::get_bits(x->entry_type_);
//     }
//
//
//
//     static et_node_ptr get_tree_edge_entry(const node_ptr& x) {
//       return compression::get_pointer<et_node_ptr>(x->entry_type_);      
//     }
//
//     static void set_tree_edge_entry(const node_ptr& x, const et_node_ptr& y) {
//       compression::set_pointer_and_bits(x->entry_type_, y, true);      
//     }
//
//
//
//     static etnte_node_ptr get_non_tree_edge_entry(const node_ptr& x) {
//       return compression::get_pointer<etnte_node_ptr>(x->entry_type_);      
//     }
//
//     static void set_non_tree_edge_entry(const node_ptr& x, const etnte_node_ptr& y) {
//       compression::set_pointer_and_bits(x->entry_type_, y, false);      
//     }    
//   };
// }