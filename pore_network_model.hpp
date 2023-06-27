#pragma once


#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/static_vector.hpp>
#include <dpl/soa.hpp>

#include <boost/iostreams/device/mapped_file.hpp>

#include <concepts>
#include <filesystem>
// #include <fstream>
#include <vector>
#include <cstring>

#include <stdint.h>






namespace xpm
{
  using v3i = dpl::vector3i;
  using v3d = dpl::vector3d;

  namespace voxel_tag
  {
    struct phase
    {
      std::uint8_t value;

      friend bool operator==(const phase& lhs, const phase& rhs) { return lhs.value == rhs.value; }
      friend bool operator!=(const phase& lhs, const phase& rhs) { return !(lhs == rhs); }
    };

    /**
     * \brief
     *   -2: non-pore (solid or microporous) | the same
     *   -1: inlet/outlet (do not know?) | the same
     *   >=0: pore clusters
     */
    struct velem
    {
      std::int32_t value;

      bool non_pore() {
        return value == -2;
      }
    };


  }

  namespace presets
  {
    static inline constexpr voxel_tag::phase pore = {0};
    static inline constexpr voxel_tag::phase solid = {1};
    static inline constexpr voxel_tag::phase microporous = {2};  
  }

  namespace geometric_properties
  {
    struct equilateral_triangle_properties
    {
      static constexpr double area(double r_ins = 1) {
        return 5.19615242271*r_ins*r_ins;
      }

      // k * G, k - coefficient, G - shape factor
      static constexpr double conductance(double area = 1, double viscosity = 1) {
        return 0.0288675134595*area*area/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }
    };
  }
  
  namespace attribs {
    def_static_key(pos)
    def_static_key(r_ins)
    def_static_key(adj)
    def_static_key(length)
    def_static_key(length0)
    def_static_key(length1)
  }

  using mapped_file_source = boost::iostreams::mapped_file_source;
  
  using pnm_idx = int64_t;
  using pnm_3idx = dpl::vector_n<pnm_idx, 3>;
  
  class pore_network_model
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
    
    
    pnm_idx parse_statoil_text_idx(std::integral auto idx) const {
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
    pnm_idx node_count_ = 0;
    pnm_idx throat_count_ = 0;


    dpl::soa<
      attribs::pos_t, dpl::vector3d,
      attribs::r_ins_t, double
    > node_;

    dpl::soa<
      attribs::adj_t, std::pair<pnm_idx, pnm_idx>,
      attribs::r_ins_t, double,
      attribs::length_t, double,
      attribs::length0_t, double,
      attribs::length1_t, double
    > throat_;
    
    

    
    // std::vector<std::pair<pnm_idx, pnm_idx>> throats;        

    // std::vector<dpl::vector3d> node_pos;

    
    // Throat values go first, then node values.      
    // std::vector<double> r_ins_; 
    // std::vector<double> volume;
    // std::vector<double> shapeFactor;
    // std::vector<double> lengthThroat;
    // std::vector<double> _length0;
    // std::vector<double> _length1;

    enum class file_format
    {
      statoil,
      binary_dp69
    };

    pore_network_model() = default;
    
    pore_network_model(const pore_network_model& other) = delete;
    pore_network_model(pore_network_model&& other) noexcept = default;
    pore_network_model& operator=(const pore_network_model& other) = delete;
    pore_network_model& operator=(pore_network_model&& other) noexcept = default;

    pore_network_model(const std::filesystem::path& p, file_format ff) {
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
    bool inner_node(pnm_idx i) const {
      return i < node_count_;
    }

    auto inlet() const {
      return node_count_;
    }

    auto outlet() const {
      return node_count_ + 1;
    }
    
    // bool inlet(pnm_idx i) const {
    //   return i == inlet();
    // }
    //
    // bool outlet(pnm_idx i) const {
    //   return i == outlet();
    // }

    // auto node_r_ins(pnm_idx i) const {
    //   return node_[attribs::r_ins][i];
    //   // return r_ins_[throat_count_ + i];
    // }

    void read_from_binary_file(const std::filesystem::path& network_path) {
      using namespace attribs;

      mapped_file_source file(network_path.string());
      auto* ptr = const_cast<char*>(file.data());

      parse_bin(ptr, node_count_);
      parse_bin(ptr, throat_count_);

      node_.resize(node_count_);
      throat_.resize(throat_count_);

      parse_bin(ptr, node_.ptr(pos), node_count_);
      parse_bin(ptr, node_.ptr(r_ins), node_count_);
      
      for (pnm_idx i = 0; i < throat_count_; ++i) {
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
        mapped_file_source node1_file(network_path_str + "_node1.dat");
        auto* node1_ptr = const_cast<char*>(node1_file.data());
        parse_text(node1_ptr, node_count_);
        parse_text(node1_ptr, physical_size);

        node_.resize(node_count_);
        
        for (pnm_idx i = 0, count = node_count_; i < count; ++i) {
          skip_word(node1_ptr);
          parse_text(node1_ptr, node_[pos][i]);
          skip_line(node1_ptr);
        }
      }

      {
        mapped_file_source throat1_file(network_path_str + "_link1.dat");
        auto* throat1_ptr = const_cast<char*>(throat1_file.data());
        parse_text(throat1_ptr, throat_count_);   

        throat_.resize(throat_count_);
        
        pnm_idx value;
        
        for (pnm_idx i = 0; i < throat_count_; ++i) {       
          skip_word(throat1_ptr);

          parse_text(throat1_ptr, value);
          throat_[adj][i].first = parse_statoil_text_idx(value); 

          parse_text(throat1_ptr, value);
          throat_[adj][i].second = parse_statoil_text_idx(value);        

          parse_text(throat1_ptr, throat_[r_ins][i]);
          // parse_text(throat1_ptr, shapeFactor[i]); //TODO

          skip_line(throat1_ptr);
        }        
      }      

      {
        mapped_file_source throat2_file(network_path_str + "_link2.dat");
        auto* throat2_ptr = const_cast<char*>(throat2_file.data());

        for (pnm_idx i = 0; i < throat_count_; ++i) {
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
        mapped_file_source node2_file(network_path_str + "_node2.dat");
        auto* node2_ptr = const_cast<char*>(node2_file.data());        

        // auto r_ins_node = r_ins_.begin() + throat_count_;
        // auto volume_node = volume.begin() + throat_count_;
        // auto shape_factor_node = shapeFactor.begin() + throat_count_;

        for (pnm_idx i = 0, count = node_count_; i < count; i++) {          
          skip_word(node2_ptr);
          skip_word(node2_ptr); // parse_text(node2_ptr, volume_node[i]); // TODO
          parse_text(node2_ptr, node_[r_ins][i]);
          // parse_text(node2_ptr, shape_factor_node[i]); // TODO
          skip_line(node2_ptr);
        }        
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
     * \brief
     * input file value description
     *   -2: solid (validated),
     *   -1: inlet/outlet (do not know?),
     *   0, 1: do not exist (validated),
     *   >=2: cluster (0, 1 are inlet and outlet node indices in Statoil format)
     */
    static std::unique_ptr<voxel_tag::velem[]> read_icl_velems(const std::filesystem::path& network_path, const pnm_3idx& dim) {
      mapped_file_source file(network_path.string() + "_VElems.raw");
      const auto* file_ptr = reinterpret_cast<const std::int32_t*>(file.data());

      auto result = std::make_unique<voxel_tag::velem[]>(dim.prod());
      auto* ptr = result.get();

      pnm_3idx velems_factor{1, dim.x() + 2, (dim.x() + 2)*(dim.y() + 2)};
      pnm_3idx ijk;
      
      for (ijk.z() = 0; ijk.z() < dim.z(); ++ijk.z())
        for (ijk.y() = 0; ijk.y() < dim.y(); ++ijk.y())
          for (ijk.x() = 0; ijk.x() < dim.x(); ++ijk.x()) {
            auto val = file_ptr[velems_factor.dot(ijk + 1)];

            *ptr++ = {val > 0 ? val - 2 : val};
          }

      return result;
    }





    struct parallel_parts_mapping
    {
      std::unique_ptr<pnm_idx[]> normal_to_optimised;
      std::unique_ptr<pnm_idx[]> optimised_to_normal;
      std::vector<std::pair<HYPRE_BigInt, HYPRE_BigInt>> rows_per_block;
    };
    
    parallel_parts_mapping parallel_partitioning(v3i blocks) {
      auto block_size = physical_size/blocks;

      v3i map_idx{1, blocks.x(), blocks.x()*blocks.y()};

      using idx_block = std::pair<pnm_idx, int>;

      std::vector<idx_block> partioning;

      for (pnm_idx i = 0; i < node_count_; ++i) {
        auto block_idx = node_[attribs::pos][i]/block_size;

        partioning.emplace_back(
          i,
          map_idx.dot({
            std::clamp<int>(std::floor(block_idx.x()), 0, blocks.x() - 1),    // NOLINT(clang-diagnostic-float-conversion)
            std::clamp<int>(std::floor(block_idx.y()), 0, blocks.y() - 1),    // NOLINT(clang-diagnostic-float-conversion)
            std::clamp<int>(std::floor(block_idx.z()), 0, blocks.z() - 1)})   // NOLINT(clang-diagnostic-float-conversion)
        );  
      }

      std::ranges::sort(partioning, [](const idx_block& l, const idx_block& r) {
        return l.second < r.second;
      });


      parallel_parts_mapping mapping;

      mapping.normal_to_optimised = std::make_unique<pnm_idx[]>(node_count_);
      mapping.optimised_to_normal = std::make_unique<pnm_idx[]>(node_count_);

      // std::vector<pnm_idx> mapping;

      int block_count = blocks.prod();
      auto row_count_per_block = std::make_unique<pnm_idx[]>(block_count);
      mapping.rows_per_block.resize(block_count);

      for (auto i = 0; i < block_count; ++i)
        row_count_per_block[i] = 0;

      for (pnm_idx i = 0; i < node_count_; ++i) {
        mapping.normal_to_optimised[partioning[i].first] = i;
        mapping.optimised_to_normal[i] = partioning[i].first;
        ++row_count_per_block[partioning[i].second];
      }

      int first_row = 0;
      for (auto i = 0; i < block_count; ++i) {
        mapping.rows_per_block[i] = {first_row, first_row + row_count_per_block[i]};
        first_row += row_count_per_block[i] + 1;
      }

      return mapping;
    }

    auto generate_pressure_input(const parallel_parts_mapping& mapping) {

      using eq_tri = geometric_properties::equilateral_triangle_properties;
      using namespace attribs;

      auto* map = mapping.normal_to_optimised.get();
      
      dpl::hypre::ls_known_storage lks;

      dpl::hypre::ls_known_storage_builder builder{lks};

      builder.allocate_rows(node_count_);

      // lks.nrows = node_count_;
      // lks.ncols = std::make_unique<HYPRE_Int[]>(node_count_);
      // lks.b = std::make_unique<HYPRE_Complex[]>(node_count_);
      //
      // for (pnm_idx i = 0; i < node_count_; ++i) {
      //   lks.ncols[i] = 1; // diag coef
      //   lks.b[i] = 0;
      // }
      //
      // lks.nvalues = node_count_; // diag coef

      for (pnm_idx i = 0; i < throat_count_; ++i) {
        auto [n0, n1] = throat_[adj][i];

        if (inner_node(n0) && inner_node(n1)) {
          builder.reserve_connection(map[n0], map[n1]);
          // ++lks.ncols[map[n0]];
          // ++lks.ncols[map[n1]];
          // lks.nvalues += 2;
        }
      }

      builder.allocate_values();

      // lks.cols = std::make_unique<HYPRE_BigInt[]>(lks.nvalues);
      // lks.values = std::make_unique<HYPRE_Complex[]>(lks.nvalues);
      //
      //
      // auto diag_shift = std::make_unique<pnm_idx[]>(node_count_);
      // auto off_shift = std::make_unique<pnm_idx[]>(node_count_);
      //
      //
      // diag_shift[0] = 0;
      // off_shift[0] = 0; // pre increment will be required
      // lks.cols[0] = 0;
      // lks.values[0] = 0; // diag
      //
      // for (pnm_idx i = 1; i < node_count_; ++i) {
      //   auto s = diag_shift[i - 1] + lks.ncols[i - 1];
      //   diag_shift[i] = s;
      //   off_shift[i] = s; // pre increment will be required
      //   lks.cols[s] = i;
      //   lks.values[s] = 0; // diag
      // }



      for (pnm_idx i = 0; i < throat_count_; ++i) {
        auto [n0, n1] = throat_[adj][i];

        // if (i%(throat_count_/10) == 0) {
        //   std::cout << (1.0*i)/throat_count_ << '\n';
        // }
        
        auto coef = -1.0/(
          (inner_node(n0) ? throat_[length0][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][n0])) : 0.0) +
          (inner_node(n1) ? throat_[length1][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][n1])) : 0.0) +
          throat_[length][i]/eq_tri::conductance(eq_tri::area(throat_[r_ins][i])));

        if (n0 == inlet()) {
          builder.add_b(map[n1], coef/**1 Pa*/);
          builder.add_diag(map[n1], coef);

          // lks.b[map[n1]] += coef/**1 Pa*/;
          // lks.values[diag_shift[map[n1]]] += coef; // diag
        }
        else if (n1 == inlet()) {
          builder.add_b(map[n0], coef/**1 Pa*/);
          builder.add_diag(map[n0], coef);

          // lks.b[map[n0]] += coef/**1 Pa*/;
          // lks.values[diag_shift[map[n0]]] += coef; // diag
        }
        else if (n0 == outlet()) {
          // builder.add_b(map[n1], coef/**0 Pa*/);
          builder.add_diag(map[n1], coef);

          // // free_terms[n1] += coef/**0 Pa*/;
          // lks.values[diag_shift[map[n1]]] += coef; // diag
        }
        else if (n1 == outlet()) {
          // builder.add_b(map[n0], coef/**0 Pa*/);
          builder.add_diag(map[n0], coef);

          // // free_terms[n0] += coef/**0 Pa*/;
          // lks.values[diag_shift[map[n0]]] += coef; // diag
        }
        else /*if (n0 != outlet() && n1 != outlet())*/ {
          builder.set_connection(map[n0], map[n1], coef);

          // lks.values[diag_shift[map[n0]]] += coef; // diag
          // auto off_n0 = ++off_shift[map[n0]];
          // lks.cols[off_n0] = map[n1];
          // lks.values[off_n0] = -coef;
          //
          //
          // lks.values[diag_shift[map[n1]]] += coef; // diag
          // auto off_n1 = ++off_shift[map[n1]];
          // lks.cols[off_n1] = map[n0];
          // lks.values[off_n1] = -coef;
        }
      }

      return lks;










      // dpl::hypre::sparse_matrix_builder matrix(node_count_);
      //
      // std::vector<double> free_terms(node_count_, 0);
      //
      //
      //
      //
      // std::cout << "\n\nSparseMatrix creation START...\n";
      //
      // for (pnm_idx i = 0; i < throat_count_; ++i) {
      //   auto [n0, n1] = throat_[adj][i];
      //
      //   if (i%(throat_count_/10) == 0) {
      //     std::cout << (1.0*i)/throat_count_ << '\n';
      //   }
      //   
      //   auto coef = -1.0/(
      //     (inner_node(n0) ? throat_[length0][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][n0])) : 0.0) +
      //     (inner_node(n1) ? throat_[length1][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][n1])) : 0.0) +
      //     throat_[length][i]/eq_tri::conductance(eq_tri::area(throat_[r_ins][i])));
      //
      //   if (n0 == inlet()) {
      //     free_terms[map[n1]] += coef/**1 Pa*/;
      //     matrix.add_diag(map[n1], coef);
      //   }
      //   else if (n1 == inlet()) {
      //     free_terms[map[n0]] += coef/**1 Pa*/;
      //     matrix.add_diag(map[n0], coef);
      //   }
      //   else if (n0 == outlet()) {
      //     // free_terms[n1] += coef/**0 Pa*/;
      //     matrix.add_diag(map[n1], coef);
      //   }
      //   else if (n1 == outlet()) {
      //     // free_terms[n0] += coef/**0 Pa*/;
      //     matrix.add_diag(map[n0], coef);
      //   }
      //   else /*if (n0 != outlet() && n1 != outlet())*/ {
      //     matrix.add_paired_diff_coef(map[n0], map[n1], -coef);
      //   }
      // }
      //
      // std::cout << "\n\nSparseMatrix creation END|||";
      //
      // // getchar();
      //
      // std::cout << "\n\nInputDeprec creation START...";
      //
      // auto inputDepr = dpl::hypre::InputDeprec /*input*/{matrix, std::move(free_terms)};
      //
      // std::cout << "\n\nInputDeprec creation END|||";
      //
      // return inputDepr;

      // return input.Solve();
    }
  };
}
