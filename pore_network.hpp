#pragma once


#include "declarations.hpp"

#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/soa.hpp>

namespace xpm
{
  class pore_network
  {
    static void skip_until(char* &ptr, char val) {   
      while (*ptr != '\0' && *ptr++ != val) {}    
    }

    static void skip_line(char* &ptr) {   
      skip_until(ptr, '\n');
    }
    
    static void skip_word(char* &ptr) {       
      while (isspace(*ptr))
        ++ptr;
      while (*ptr != '\0' && !isspace(*ptr))
        ++ptr;
    }


    // template<typename T>
    // static T parse_text(char* &ptr);
    //
    // template<>
    // static int parse_text(char* &ptr) {
    //   return strtol(ptr, &ptr, 10);
    // }
    
    static void parse_text(char* &ptr, int& val) {
      val = strtol(ptr, &ptr, 10);
    }
    
    static void parse_text(char* &ptr, long long& val) {
      val = strtoll(ptr, &ptr, 10);
    }

    static void parse_text(char* &ptr, unsigned long long& val) {
      val = strtoull(ptr, &ptr, 10);
    }

    static void parse_text(char* &ptr, double& val) {
      val = strtod(ptr, &ptr);
    }

    static void parse_text(char* &ptr, dpl::vector3d& val) {
      parse_text(ptr, val.x());
      parse_text(ptr, val.y());
      parse_text(ptr, val.z());
    }

    static void parse_bin(char* &src, auto& dst) {
      std::memcpy(&dst, src, sizeof dst);
      src += sizeof dst;
    }

    // template <typename T>
    // static void parse_bin(char* &ptr, std::vector<T>& vec, size_t count) {
    //   vec.resize(count);
    //   auto total_size = sizeof(T)*count;
    //   std::memcpy(vec.data(), ptr, total_size);
    //   ptr += total_size;
    // }

    template <typename T>
    static void parse_bin(char* &src, T* dst, size_t count) {
      auto total_size = sizeof(T)*count;
      std::memcpy(dst, src, total_size);
      src += total_size;
    }
    
    //
    // void reset() {      // TODO: Not good function
    //   throats.resize(throat_count_);
    //   lengthThroat.resize(throat_count_);
    //   _length0.resize(throat_count_, 0);
    //   _length1.resize(throat_count_, 0);
    //
    //
    //   auto inlet_total_idx = throat_count_ + node_count_ - 2;
    //   auto outlet_total_idx = inlet_total_idx + 1;
    //   r_ins_.resize(throat_count_ + node_count_);            
    //   r_ins_[inlet_total_idx] = 1e30;
    //   r_ins_[outlet_total_idx] = 1e-30;
    //
    //   volume.resize(throat_count_ + node_count_);
    //   volume[inlet_total_idx] = 0.0;
    //   volume[outlet_total_idx] = 0.0;
    //
    //   shapeFactor.resize(throat_count_ + node_count_);
    //   shapeFactor[inlet_total_idx] = -9999.0;
    //   shapeFactor[outlet_total_idx] = -9999.0;
    // }

    //
    // void default_inlet_outlet() {      // TODO: Not good function
    //   auto inlet_total_idx = node_count_ - 2;
    //   auto outlet_total_idx = inlet_total_idx + 1;
    //   
    //   node_[attribs::r_ins][inlet_total_idx] = 1e30;
    //   node_[attribs::r_ins][outlet_total_idx] = 1e-30;
    //
    //   // volume.resize(throat_count_ + node_count_);  // TODO
    //   // volume[inlet_total_idx] = 0.0;
    //   // volume[outlet_total_idx] = 0.0;
    //
    //   // shapeFactor.resize(throat_count_ + node_count_);
    //   // shapeFactor[inlet_total_idx] = -9999.0;
    //   // shapeFactor[outlet_total_idx] = -9999.0;
    // }
    
    
    idx1d_t parse_statoil_text_idx(std::integral auto idx) const {
      if (idx == -1)
        return node_count_;
  
      if (idx == 0)
        return node_count_ + 1;
    
      return idx - 1;
    }

    
  
  public:
    dpl::vector3d physical_size{1};

    /**
     * \brief
     *    inner node count \n
     *    disregards inlet and outlet \n
     *    INLET_IDX = node_count_ \n
     *    OUTLET_IDX = node_count_ + 1
     */
    idx1d_t node_count_ = 0;
    size_t throat_count_ = 0;


    dpl::soa<
      attribs::pos_t, dpl::vector3d,
      attribs::r_ins_t, double
    > node_;

    /**
     * \brief
     * (left < right) in adj_t\n
     * left is always an inner node\n
     * right is an inner or outer (inlet/outlet) node
     */
    dpl::soa<
      attribs::adj_t, std::pair<idx1d_t, idx1d_t>,     
      attribs::r_ins_t, double,
      attribs::length_t, double,
      attribs::length0_t, double,
      attribs::length1_t, double
    > throat_;
    

    enum class file_format
    {
      statoil,
      binary_dp69
    };

    pore_network() = default;
    ~pore_network() = default;

    pore_network(const pore_network& other) = delete;
    pore_network(pore_network&& other) noexcept = default;
    pore_network& operator=(const pore_network& other) = delete;
    pore_network& operator=(pore_network&& other) noexcept = default;

    pore_network(const std::filesystem::path& p, file_format ff) {
      if (ff == file_format::statoil)
        read_from_text_file(p);
      else if (ff == file_format::binary_dp69)
        read_from_binary_file(p);
    }

    


    
    
    // auto node_count() const {
    //   return node_count_;
    // }

    // auto throat_count() const {
    //   return throat_count_;
    // }

    // auto inner_count() const {
    //   return inner_node_count() + throat_count();
    // }

    void scale(double x) {
      using namespace attribs;
            
      for (auto& p : node_.range(pos))
        p = p*x;

      for (auto& r : node_.range(r_ins))
        r = r*x;

      for (auto& r : throat_.range(r_ins))
        r = r*x;
    }
    
    /**
     * \brief inlet and outlet are outer nodes, other nodes are inner.
     */
    bool inner_node(idx1d_t i) const {
      return i < node_count_;
    }

    auto inlet() const {
      return node_count_;
    }

    auto outlet() const {
      return node_count_ + 1;
    }
    
    // auto node_r_ins(pnm_idx i) const {
    //   return node_[attribs::r_ins][i];
    //   // return r_ins_[throat_count_ + i];
    // }

    void read_from_binary_file(const std::filesystem::path& network_path) { // TODO
      throw std::exception("[read_from_binary_file] not implemented correctly");

      using namespace attribs;

      boost::iostreams::mapped_file_source file(network_path.string());
      auto* ptr = const_cast<char*>(file.data());

      parse_bin(ptr, node_count_);
      parse_bin(ptr, throat_count_); // TODO

      node_.resize(node_count_);
      throat_.resize(throat_count_);

      parse_bin(ptr, node_.ptr(pos), node_count_);
      parse_bin(ptr, node_.ptr(r_ins), node_count_);
      
      for (idx1d_t i = 0; i < throat_count_; ++i) {
        dpl::vector2i pair;
        parse_bin(ptr, pair);
        throat_[adj][i] = {pair.x(), pair.y()};
      }

      parse_bin(ptr, throat_.ptr(r_ins), throat_count_);

      for (auto& r : node_.range(r_ins))
        r = r/2;

      for (auto& r : throat_.range(r_ins))
        r = r/2;
    }
    
    void read_from_text_file(const std::filesystem::path& network_path) {
      using namespace attribs;

      auto network_path_str = network_path.string();

      {
        boost::iostreams::mapped_file_source node1_file(network_path_str + "_node1.dat");
        auto* node1_ptr = const_cast<char*>(node1_file.data());
        parse_text(node1_ptr, node_count_);
        parse_text(node1_ptr, physical_size);

        node_.resize(node_count_);
        
        for (idx1d_t i = 0, count = node_count_; i < count; ++i) {
          skip_word(node1_ptr);
          parse_text(node1_ptr, node_[pos][i]);
          skip_line(node1_ptr);
        }
      }

      {
        boost::iostreams::mapped_file_source throat1_file(network_path_str + "_link1.dat");
        auto* throat1_ptr = const_cast<char*>(throat1_file.data());
        parse_text(throat1_ptr, throat_count_);   

        throat_.resize(throat_count_);
        
        idx1d_t value;
        
        for (size_t i = 0; i < throat_count_; ++i) {       
          skip_word(throat1_ptr);

          auto& [left, right] = throat_[adj][i];

          parse_text(throat1_ptr, value);
          left = parse_statoil_text_idx(value); 

          parse_text(throat1_ptr, value);
          right = parse_statoil_text_idx(value);        
          
          parse_text(throat1_ptr, throat_[r_ins][i]);
          // parse_text(throat1_ptr, shapeFactor[i]); //TODO

          skip_line(throat1_ptr);
        }        
      }      

      {
        boost::iostreams::mapped_file_source throat2_file(network_path_str + "_link2.dat");
        auto* throat2_ptr = const_cast<char*>(throat2_file.data());

        for (size_t i = 0; i < throat_count_; ++i) {
          dpl::sfor<3>([&throat2_ptr] {
            skip_word(throat2_ptr);  
          });
          
          parse_text(throat2_ptr, throat_[length0][i]); 
          parse_text(throat2_ptr, throat_[length1][i]);
          parse_text(throat2_ptr, throat_[length][i]);
          // parse_text(throat2_ptr, volume[i]);          //TODO
          skip_line(throat2_ptr);
        }  
      }

      {
        boost::iostreams::mapped_file_source node2_file(network_path_str + "_node2.dat");
        auto* node2_ptr = const_cast<char*>(node2_file.data());        

        // auto r_ins_node = r_ins_.begin() + throat_count_;
        // auto volume_node = volume.begin() + throat_count_;
        // auto shape_factor_node = shapeFactor.begin() + throat_count_;

        for (idx1d_t i = 0, count = node_count_; i < count; i++) {          
          skip_word(node2_ptr);
          skip_word(node2_ptr); // parse_text(node2_ptr, volume_node[i]); // TODO
          parse_text(node2_ptr, node_[r_ins][i]);
          // parse_text(node2_ptr, shape_factor_node[i]); // TODO
          skip_line(node2_ptr);
        }        
      }

      for (size_t i = 0; i < throat_count_; i++)
        if (auto& [l, r] = throat_[adj][i];
          l > r) {
          std::swap(l, r);
          std::swap(throat_[length0][i], throat_[length1][i]);
        }
    }

    // auto read_icl_velems(const std::filesystem::path& network_path) const {
    //   mapped_file_source file(network_path.string() + "_VElems.raw");
    //   auto* ptr = const_cast<char*>(file.data());
    //   auto size = file.size()/4;
    //   std::vector<std::int32_t> map(size);
    //   parse_bin(ptr, map.data(), size);
    //   return map;
    // }



    /**
     * \brief coefficient of a throat
     * \param i throat index
     */
    double coef(size_t i) {
      using eq_tri = geometric_properties::equilateral_triangle_properties;
      using namespace attribs;

      auto [left, right] = throat_[adj][i];

      return -1.0/(
        throat_[length0][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][left])) +
        throat_[length][i]/eq_tri::conductance(eq_tri::area(throat_[r_ins][i])) +
        (inner_node(right) ? throat_[length1][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][right])) : 0.0));
    }


    





    struct rows_decomposition
    {
      std::unique_ptr<idx1d_t[]> normal_to_decomposed;
      std::unique_ptr<idx1d_t[]> decomposed_to_normal;
      std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>> rows_per_block; // [from, to] inclusive in both sides
    };
    
    rows_decomposition decompose_rows(v3i blocks) {
      auto block_size = physical_size/blocks;

      using idx_block = std::pair<idx1d_t, int>;

      std::vector<idx_block> block_ordered_indices(node_count_);

      auto map_idx = idx_mapper(blocks);

      for (auto i : dpl::range(node_count_)) {
        v3i block_idx = node_[attribs::pos][i]/block_size;

        block_ordered_indices[i] = {i,
          map_idx(
            std::clamp(block_idx.x(), 0, blocks.x() - 1), // NOLINT(clang-diagnostic-float-conversion)
            std::clamp(block_idx.y(), 0, blocks.y() - 1), // NOLINT(clang-diagnostic-float-conversion)
            std::clamp(block_idx.z(), 0, blocks.z() - 1)) // NOLINT(clang-diagnostic-float-conversion)
        };
      }

      std::ranges::sort(block_ordered_indices, [](const idx_block& l, const idx_block& r) { return l.second < r.second; });

      rows_decomposition mapping;

      mapping.normal_to_decomposed = std::make_unique<idx1d_t[]>(node_count_);
      mapping.decomposed_to_normal = std::make_unique<idx1d_t[]>(node_count_);

      auto block_count = blocks.prod();
      auto row_count_per_block = std::make_unique<idx1d_t[]>(block_count);
      mapping.rows_per_block.resize(block_count);

      for (auto i : dpl::range(block_count))
        row_count_per_block[i] = 0;
      
      for (auto i : dpl::range(node_count_)) {
        mapping.normal_to_decomposed[block_ordered_indices[i].first] = i;
        mapping.decomposed_to_normal[i] = block_ordered_indices[i].first;
        ++row_count_per_block[block_ordered_indices[i].second];
      }

      int first_row = 0;
      for (auto i : dpl::range(block_count)) {
        mapping.rows_per_block[i] = {first_row, first_row + row_count_per_block[i] - 1};
        first_row += row_count_per_block[i];
      }

      return mapping;
    }


    auto generate_pressure_input_BASIC() {
      using namespace attribs;

      dpl::hypre::ls_known_storage_builder builder;

      builder.allocate_rows(node_count_);

      for (auto i : dpl::range(throat_count_))
        if (auto [l, r] = throat_[adj][i]; inner_node(r))
          builder.reserve_connection(l, r);

      builder.allocate_values();

      for (auto i : dpl::range(throat_count_)) {
        auto [l, r] = throat_[adj][i];

        auto coef = this->coef(i);

        if (r == inlet()) {
          builder.add_b(l, coef/**1 Pa*/);
          builder.add_diag(l, coef);
        }
        else if (r == outlet()) {
          // builder.add_b(map[n0], coef/**0 Pa*/);
          builder.add_diag(l, coef);
        }
        else
          builder.set_connection(l, r, coef);
      }

      return builder.acquire_storage();
    }


    auto generate_pressure_input(const rows_decomposition& mapping, const auto& coef_map) {
      using namespace attribs;

      auto* map = mapping.normal_to_decomposed.get();
      
      dpl::hypre::ls_known_storage_builder builder;

      builder.allocate_rows(node_count_);

      for (size_t i = 0; i < throat_count_; ++i)
        if (auto [n0, n1] = throat_[adj][i]; inner_node(n0) && inner_node(n1))
          builder.reserve_connection(map[n0], map[n1]);

      builder.allocate_values();

      for (size_t i = 0; i < throat_count_; ++i) {
        auto [l, r] = throat_[adj][i];

        auto coef = coef_map(i);

        if (r == inlet()) {
          builder.add_b(map[l], coef/**1 Pa*/);
          builder.add_diag(map[l], coef);
        }
        else if (r == outlet()) {
          // builder.add_b(map[n0], coef/**0 Pa*/);
          builder.add_diag(map[l], coef);
        }
        else
          builder.set_connection(map[l], map[r], coef);
      }

      return builder.acquire_storage();
    }
  };
}
