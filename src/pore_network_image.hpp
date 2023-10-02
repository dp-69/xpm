#pragma once


#include "declarations.hpp"

#include <dpl/soa.hpp>
#include <dpl/graph/dc_graph.hpp>
#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/hypre/mpi_module.hpp>

#include <boost/format.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <fmt/format.h>




namespace xpm
{
  struct rows_mapping
  {
    static inline constexpr auto invalid_block = std::numeric_limits<int>::max();

    dpl::strong_vector<net_tag, idx1d_t> forward;
    std::unique_ptr<net_idx_t[]> backward;
    std::vector<dpl::hypre::index_range> block_rows;

    rows_mapping(idx1d_t rows, int blocks) {
      forward.resize(rows);
      backward = std::make_unique<net_idx_t[]>(rows);
      block_rows = std::vector<dpl::hypre::index_range>(blocks);
    }
  };

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

    #ifdef __linux__
    static void parse_text(char* &ptr, size_t& val) { // for Ubuntu
      val = strtoull(ptr, &ptr, 10);
    }
    #endif

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
    
    
    macro_idx_t parse_statoil_text_idx(std::integral auto idx) const {
      if (idx == -1)
        return node_count();
  
      if (idx == 0)
        return node_count() + 1;
    
      return macro_idx_t{idx - 1};
    }

  public:
    dpl::vector3d physical_size{1};

    dpl::soa<
      attribs::pos_t, dpl::vector3d,
      attribs::r_ins_t, double,
      attribs::volume_t, double
    > node_;

    /**
     * \brief
     * (left < right) in adj_t\n
     * left is always an inner node\n
     * right is an inner or outer (inlet/outlet) node
     */
    dpl::soa<
      attribs::adj_t, std::pair<macro_idx_t, macro_idx_t>,     
      attribs::r_ins_t, double,
      attribs::length_t, double,
      attribs::length0_t, double,
      attribs::length1_t, double,
      attribs::volume_t, double
    > throat_;
    

    enum class file_format
    {
      statoil,
      // binary_dp69
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
      // else if (ff == file_format::binary_dp69)
      //   read_from_binary_file(p);
    }

    


    
    
    macro_idx_t node_count() const {
      return macro_idx_t(node_.size());
    }

    std::size_t throat_count() const {
      return throat_.size();
    }

    bool inner_node(macro_idx_t i) const {
      return i < node_count();
    }

    auto inlet() const {
      return macro_idx_t{node_count()};
    }

    auto outlet() const {
      return macro_idx_t{node_count() + 1};
    }

    

    


    // void read_from_binary_file(const std::filesystem::path& network_path) {
    //   using namespace attribs;
    //
    //   boost::iostreams::mapped_file_source file(network_path.string());
    //   auto* ptr = const_cast<char*>(file.data());
    //
    //   idx1d_t node_count;
    //   idx1d_t throat_count;
    //
    //   parse_bin(ptr, node_count);
    //   parse_bin(ptr, throat_count);
    //
    //   node_.resize(node_count);
    //   throat_.resize(throat_count);
    //
    //   parse_bin(ptr, node_.ptr(pos), node_count);
    //   parse_bin(ptr, node_.ptr(r_ins), node_count);
    //   
    //   for (size_t i = 0; i < throat_count; ++i) {
    //     dpl::vector2i pair;
    //     parse_bin(ptr, pair);
    //     throat_[adj][i] = {{pair.x()}, {pair.y()}};
    //   }
    //
    //   parse_bin(ptr, throat_.ptr(r_ins), throat_count);
    //
    //   for (auto& r : node_.range(r_ins))
    //     r = r/2;
    //
    //   for (auto& r : throat_.range(r_ins))
    //     r = r/2;
    // }
    
    void read_from_text_file(const std::filesystem::path& network_path) {
      using namespace attribs;

      auto network_path_str = network_path.string();

      idx1d_t node_count{0};
      size_t throat_count{0};

      {
        boost::iostreams::mapped_file_source node1_file(network_path_str + "_node1.dat");
        auto* node1_ptr = const_cast<char*>(node1_file.data());

        parse_text(node1_ptr, node_count);
        parse_text(node1_ptr, physical_size);

        node_.resize(node_count);
        
        for (idx1d_t i = 0; i < node_count; ++i) {
          skip_word(node1_ptr);
          parse_text(node1_ptr, node_[pos][i]);
          skip_line(node1_ptr);
        }
      }

      {
        boost::iostreams::mapped_file_source throat1_file(network_path_str + "_link1.dat");
        auto* throat1_ptr = const_cast<char*>(throat1_file.data());
        parse_text(throat1_ptr, throat_count);   

        throat_.resize(throat_count);
        
        idx1d_t value;
        
        for (size_t i = 0; i < throat_count; ++i) {       
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

        for (size_t i = 0; i < throat_count; ++i) {
          dpl::sfor<3>([&throat2_ptr] {
            skip_word(throat2_ptr);  
          });
          
          parse_text(throat2_ptr, throat_[length0][i]); 
          parse_text(throat2_ptr, throat_[length1][i]);
          parse_text(throat2_ptr, throat_[length][i]);
          parse_text(throat2_ptr, throat_[volume][i]);
          skip_line(throat2_ptr);
        }  
      }

      {
        boost::iostreams::mapped_file_source node2_file(network_path_str + "_node2.dat");
        auto* node2_ptr = const_cast<char*>(node2_file.data());        

        // auto r_ins_node = r_ins_.begin() + throat_count_;
        // auto volume_node = volume.begin() + throat_count_;
        // auto shape_factor_node = shapeFactor.begin() + throat_count_;

        for (idx1d_t i = 0; i < node_count; ++i) {          
          skip_word(node2_ptr);
          parse_text(node2_ptr, node_[volume][i]);
          parse_text(node2_ptr, node_[r_ins][i]);
          // parse_text(node2_ptr, shape_factor_node[i]); // TODO
          skip_line(node2_ptr);
        }        
      }

      for (size_t i = 0; i < throat_count; ++i)
        if (auto& [l, r] = throat_[adj][i];
          r < l) {
          std::swap(l, r);
          std::swap(throat_[length0][i], throat_[length1][i]);
        }
    }


    double coef(std::size_t i) const {
      using eq_tri = hydraulic_properties::equilateral_triangle_properties;
      using namespace attribs;

      auto [left, right] = throat_[adj][i];

      return -1.0/(
        throat_[length0][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][*left])) +
        throat_[length][i]/eq_tri::conductance(eq_tri::area(throat_[r_ins][i])) +
        (inner_node(right) ? throat_[length1][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][*right])) : 0.0));
    }

    bool eval_inlet_outlet_connectivity() const {
      disjoint_sets ds(*node_count() + 2);

      for (auto [l, r] : throat_.span(attribs::adj))
        ds.union_set(*l, *r);

      return ds.find_set(*inlet()) == ds.find_set(*outlet());
    }

    std::pair<dpl::hypre::ls_known_storage, size_t> generate_pressure_input() const {
      dpl::hypre::ls_known_storage_builder builder{*node_count()};

      for (std::size_t i = 0; i < throat_count(); ++i)
        if (auto [l, r] = throat_[attribs::adj][i]; inner_node(r))
          builder.reserve(*l, *r);

      builder.allocate();

      for (std::size_t i = 0; i < throat_count(); ++i) {
        auto [l, r] = throat_[attribs::adj][i];

        auto coef = this->coef(i);

        if (r == inlet()) {
          builder.add_b(*l, coef/**1 Pa*/);
          builder.add_diag(*l, coef);
        }
        else if (r == outlet()) {
          // builder.add_b(map[n0], coef/**0 Pa*/);
          builder.add_diag(*l, coef);
        }
        else
          builder.set(*l, *r, coef);
      }

      return {builder.acquire(), builder.nvalues()};
    }


     
    void connectivity_flow_summary(const dpl::hypre::mpi::solve_result& solve  /*const std::unique_ptr<double[]>& pressure*/) const {
      auto& [pressure, residual, iter] = solve; 

      disjoint_sets ds(*node_count());
      std::vector<bool> connected_inlet(*node_count());
      std::vector<bool> connected_outlet(*node_count());
      
      auto connected = [&](macro_idx_t i) {
        auto rep = ds.find_set(*i);
        return connected_inlet[rep] && connected_outlet[rep];
      };
      
      {
        for (auto [l, r] : throat_.span(attribs::adj))
          if (inner_node(r))
            ds.union_set(*l, *r);
      
        for (auto [l, r] : throat_.span(attribs::adj))
          if (r == inlet())
            connected_inlet[ds.find_set(*l)] = true;
          else if (r == outlet())
            connected_outlet[ds.find_set(*l)] = true;
      
        idx1d_t disconnected_macro = 0;
        
        for (macro_idx_t i{0}; i < node_count(); ++i)
          if (!connected(i))
            ++disconnected_macro;
      
        double inlet_flow = 0;
        double outlet_flow = 0;
      
        std::cout << fmt::format("  isolated {} nodes", disconnected_macro);
      
        for (std::size_t i = 0; i < throat_count(); ++i) {
          auto [l, r] = throat_[attribs::adj][i];
      
          if (r == inlet())
            if (connected(l))
              inlet_flow += coef(i)*(1 - pressure[*l]);
      
          if (r == outlet())
            if (connected(l))
              outlet_flow += coef(i)*(pressure[*l]);
        }
      
        std::cout << fmt::format(R"(
  inlet perm: {:.6f} mD
  outlet perm: {:.6f} mD
  residual: {:.4g}, iterations: {}
)",
          -inlet_flow/physical_size.x()/presets::darcy_to_m2*1000,
          -outlet_flow/physical_size.x()/presets::darcy_to_m2*1000,
          residual,
          iter);
      }
    }

    void connectivity_flow_summary_MPI(HYPRE_Real tolerance, HYPRE_Int max_iterations) const {
      auto [input, nvalues] = generate_pressure_input();
      auto nrows = *node_count();

      dpl::hypre::mpi::save(input, nrows, nvalues, {{0, nrows - 1}}, tolerance, max_iterations);

      std::system(fmt::format("mpiexec -np 1 \"{}\" -s",  // NOLINT(concurrency-mt-unsafe)
        dpl::hypre::mpi::mpi_exec).c_str());

      connectivity_flow_summary(dpl::hypre::mpi::load_values(nrows));
    }


    void connectivity_flow_summary(HYPRE_Real tolerace, HYPRE_Int max_iterations) const {
      auto pressure = std::make_unique<double[]>(*node_count());
      auto [residual, iters] = solve(
        generate_pressure_input().first, {0, *node_count() - 1}, pressure.get(), tolerace, max_iterations); // gross solve (with isolated)
      connectivity_flow_summary({std::move(pressure), residual, iters});
    }
  };




  struct image_data
  {
    idx3d_t dim;
    idx1d_t size;

    dpl::strong_vector<voxel_tag, voxel_property::phase_t> phase;

    /**
     * \brief
     *    for a void voxel - a macro node it belongs
     *    for a microporous voxel - an adjacent pore node
     */
    dpl::strong_vector<voxel_tag, voxel_property::velem_t> velem;

    auto idx1d_mapper() const {
      return idx_mapper(dim);
    }

    void eval_microporous() {
      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      voxel_idx_t idx1d{0};

      auto map_idx = idx1d_mapper();

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i, ++idx1d)
            if (phase[idx1d] == presets::microporous) {
              voxel_property::velem_t adj;

              dpl::sfor<3>([&](auto d) {
                if (ijk[d] > 0)
                  if (voxel_idx_t adj_idx = idx1d - map_idx[d]; phase[adj_idx] == presets::pore)
                    if (velem[adj_idx])
                      adj = velem[adj_idx];

                if (ijk[d] < dim[d] - 1)
                  if (voxel_idx_t adj_idx = idx1d + map_idx[d]; phase[adj_idx] == presets::pore)
                    if (velem[adj_idx])
                      adj = velem[adj_idx];
              });

              velem[idx1d] = adj;
            }
    }

    void read_image(const auto& image_path, parse::image_dict input_config) {
      std::ifstream is(image_path);
      is.seekg(0, std::ios_base::end);
      size = static_cast<idx1d_t>(is.tellg());
      is.seekg(0, std::ios_base::beg);
      phase.resize(size);
      is.read(reinterpret_cast<char*>(phase.data()), size);

      size_t pore_voxels = 0;
      size_t solid_voxels = 0;
      size_t microporous_voxels = 0;

      for (voxel_idx_t i{0}; i < size; ++i)
        if (auto& value = phase[i];
          *value == input_config.pore) {
          value = presets::pore;
          ++pore_voxels;
        }
        else if (*value == input_config.solid) {
          value = presets::solid;
          ++solid_voxels;
        }
        else if (*value == input_config.microporous) { 
          value = presets::microporous;
          ++microporous_voxels;
        }

      // #ifdef XPM_DEBUG_OUTPUT // NOLINTNEXTLINE(clang-diagnostic-misleading-indentation)
      std::cout << fmt::format(
R"(image voxels
  total: {:L}
  pore: {:L}
  solid: {:L}
  microprs: {:L}

)",
        size, pore_voxels, solid_voxels, microporous_voxels); 
      // #endif
    }

    /**
     * \brief
     * input file value description
     *   -2: solid (validated),
     *   -1: inlet/outlet (do not know?),
     *   0, 1: do not exist (validated),
     *   >=2: cluster (0, 1 are inlet and outlet node indices in Statoil format)
     */
    void read_icl_velems(const std::filesystem::path& network_path) {
      boost::iostreams::mapped_file_source file(network_path.string() + "_VElems.raw");
      const auto* file_ptr = reinterpret_cast<const std::int32_t*>(file.data());

      velem.resize(dim.prod());

      // auto* ptr = velem.get();

      idx3d_t velems_factor{1, dim.x() + 2, (dim.x() + 2)*(dim.y() + 2)};
      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      voxel_idx_t idx1d{0};

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i, ++idx1d) {
            if (auto val = file_ptr[velems_factor.dot(ijk + 1)]; val > 0)
              *velem[idx1d] = val - 2;
          }
    }
  };



  class pore_network_image
  {
    pore_network* pn_;
    image_data* img_;

    static inline constexpr net_idx_t isolated_idx_{std::numeric_limits<idx1d_t>::max()};

    struct total_tag {};
    using total_idx_t = dpl::strong_integer<idx1d_t, total_tag>;

    /**
     * \brief
     *   .size() = pn_->node_count_ + img_->size
     *   maps total() index to compressed effective (flowing) index
     */
    dpl::strong_vector<total_tag, net_tag, idx1d_t> total_to_net_map_;
    dpl::strong_vector<net_tag, total_tag, idx1d_t> net_to_total_map_;

    /**
     * \brief connected/flowing
     */
    net_idx_t connected_macro_count_{0};
    net_idx_t connected_count_{0};

    total_idx_t total(macro_idx_t i) const { return total_idx_t{*i}; }
    total_idx_t total(voxel_idx_t i) const { return total_idx_t{*pn_->node_count() + *i}; }

  public:
    pore_network_image(pore_network& pn, image_data& img)
      : pn_(&pn), img_(&img) {}

    auto& pn() const {
      return *pn_;
    }

    auto& img() const {
      return *img_;
    }

    net_idx_t connected_macro_count() const { return connected_macro_count_; }
    net_idx_t connected_count() const { return connected_count_; }

    net_idx_t net(macro_idx_t i) const { return total_to_net_map_[total(i)]; }
    net_idx_t net(voxel_idx_t i) const { return total_to_net_map_[total(i)]; }

    bool is_macro(net_idx_t i) const {
      return i < connected_macro_count_;
    }

    macro_idx_t macro(net_idx_t i) const {
      return macro_idx_t{*net_to_total_map_[i]};
    }

    voxel_idx_t voxel(net_idx_t i) const {
      return voxel_idx_t{*net_to_total_map_[i] - *pn_->node_count()};
    }

    bool connected(macro_idx_t i) const {
      return net(i) != isolated_idx_;
    }

    /**
     * \brief if connected, the voxel must be microporous, i.e. the latter evaluation is redundant
     */
    bool connected(voxel_idx_t i) const {
      return net(i) != isolated_idx_;
    }

    void evaluate_isolated() {
      auto gross_total_size = *pn_->node_count() + img_->size;

      disjoint_sets ds(gross_total_size);

      for (auto [l, r] : pn_->throat_.span(attribs::adj))
        if (pn_->inner_node(r)) // macro-macro
          ds.union_set(*l, *r);

      using namespace presets;
      auto map_idx = img_->idx1d_mapper();

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (img_->phase[idx1d] == microporous) {
                if (img_->velem[idx1d]) // macro-darcy
                  ds.union_set(*total(img_->velem[idx1d]), *total(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx_t adj_idx1d = idx1d + map_idx[d]; img_->phase[adj_idx1d] == microporous) // darcy-darcy
                      ds.union_set(*total(idx1d), *total(adj_idx1d));
                });
              }
      }

      std::vector<bool> inlet(gross_total_size);
      std::vector<bool> outlet(gross_total_size);

      for (auto [l, r] : pn_->throat_.span(attribs::adj))
        if (r == pn_->inlet()) // macro-inlet
          inlet[ds.find_set(*total(l))] = true;
        else if (r == pn_->outlet()) // macro-outlet
          outlet[ds.find_set(*total(l))] = true;

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx_t inlet_idx1d{map_idx(0, j, k)}; img_->phase[inlet_idx1d] == microporous) // darcy-inlet
              inlet[ds.find_set(*total(inlet_idx1d))] = true;

            if (voxel_idx_t outlet_idx1d{map_idx(img_->dim.x() - 1, j, k)}; img_->phase[outlet_idx1d] == microporous) // darcy-outlet
              outlet[ds.find_set(*total(outlet_idx1d))] = true;
          }
      }

      total_to_net_map_.resize(gross_total_size);

      for (macro_idx_t i{0}; i < pn_->node_count(); ++i) {
        auto total_idx = total(i);
        auto rep = ds.find_set(*total_idx); 
        total_to_net_map_[total_idx] = inlet[rep] && outlet[rep] ? connected_macro_count_++ : isolated_idx_;
      }
      
      connected_count_ = connected_macro_count_;

      for (voxel_idx_t i{0}; i < img_->size; ++i) {
        auto total_idx = total(i);
        auto rep = ds.find_set(*total_idx); 
        total_to_net_map_[total_idx] = inlet[rep] && outlet[rep] ? connected_count_++ : isolated_idx_;
      }

      net_to_total_map_.resize(*connected_count_);

      for (total_idx_t i{0}; i < gross_total_size; ++i)
        if (total_to_net_map_[i] != isolated_idx_)
          net_to_total_map_[total_to_net_map_[i]] = i;
    }


    std::tuple<idx1d_t, rows_mapping> generate_mapping(const dpl::vector3i& blocks, auto filter) const {
      auto block_size = pn_->physical_size/blocks;

      using pair = std::pair<net_idx_t, int>;

      std::vector<pair> idx_to_block(*connected_count_);

      auto map_idx = idx_mapper(blocks);

      auto eval_block = [&](const dpl::vector3d& pos) {
        dpl::vector3i block_idx = pos/block_size;

        return *map_idx(
            std::clamp(block_idx.x(), 0, blocks.x() - 1), 
            std::clamp(block_idx.y(), 0, blocks.y() - 1), 
            std::clamp(block_idx.z(), 0, blocks.z() - 1));
      };

      {
        for (macro_idx_t i{0}; i < pn_->node_count(); ++i)
          if (connected(i)) { // macro node
            net_idx_t net_idx = net(i); 
            idx_to_block[*net_idx] = {net_idx, filter(net_idx)
              ? eval_block(pn_->node_[attribs::pos][*i])
              : rows_mapping::invalid_block};
          }

        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        auto cell_size = pn_->physical_size/img_->dim;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) { // darcy node
                net_idx_t net_idx = net(idx1d);
                idx_to_block[*net_idx] = {net_idx, filter(net_idx)
                  ? eval_block(cell_size*(ijk + 0.5))
                  : rows_mapping::invalid_block};
              }
      }

      std::ranges::sort(idx_to_block, [](const pair& l, const pair& r) { return l.second < r.second; });

      auto end = std::ranges::lower_bound(idx_to_block, rows_mapping::invalid_block, {}, [](const pair& p) { return p.second; });
      idx1d_t filtered_count = end - idx_to_block.begin();
      // idx_to_block.resize(filtered_count);

      auto block_count = blocks.prod();

      rows_mapping mapping{*connected_count_, block_count};
      
      auto count_per_block = std::make_unique<idx1d_t[]>(block_count);

      for (int i = 0; i < block_count; ++i)
        count_per_block[i] = 0;
      
      for (idx1d_t i{0}; i < *connected_count_; ++i) {
        mapping.forward[idx_to_block[i].first] = i;
        mapping.backward[i] = idx_to_block[i].first;
        if (idx_to_block[i].second != rows_mapping::invalid_block)
          ++count_per_block[idx_to_block[i].second];
      }

      HYPRE_BigInt first_row = 0;
      for (int i = 0; i < block_count; ++i) {
        mapping.block_rows[i] = {first_row, first_row + count_per_block[i] - 1};
        first_row += count_per_block[i];
      }

      return {filtered_count, std::move(mapping)};
    }


    template <bool Outlet = false>
    std::tuple<
      dpl::graph::dc_graph,
      std::unordered_map<dpl::graph::dc_graph::edge_descriptor, std::size_t>,
      std::unique_ptr<dpl::graph::dc_graph::edge_descriptor[]>>
    generate_dc_graph() const {
      using namespace dpl::graph;
      using namespace presets;
      auto map_idx = img_->idx1d_mapper();

      net_idx_t vertex_count = connected_count_; // last is outlet

      if constexpr (Outlet)
        ++vertex_count;

      graph_generator<net_idx_t> gen(*vertex_count);
      std::unordered_map<dc_graph::edge_descriptor, std::size_t> de_to_throat;
      auto throat_to_de = std::make_unique<dc_graph::edge_descriptor[]>(pn_->throat_count());

      for (auto [l, r] : pn_->throat_.span(attribs::adj))
        if (connected(l))
          if (pn_->inner_node(r)) // macro-macro
            gen.reserve(net(l), net(r));
          else if constexpr (Outlet) {
            if (pn_->outlet() == r) // macro-outlet
              gen.reserve(net(l), vertex_count - 1);
          }

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if constexpr (Outlet) 
              if (voxel_idx_t adj_idx1d{map_idx(img_->dim.x() - 1, j, k)}; connected(adj_idx1d)) // darcy-outlet
                gen.reserve(net(adj_idx1d), vertex_count - 1);

            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (img_->velem[idx1d]) // macro-darcy 
                  gen.reserve(net(img_->velem[idx1d]), net(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx_t adj_idx1d = idx1d + map_idx[d]; img_->phase[adj_idx1d] == microporous) // darcy-darcy
                      gen.reserve(net(idx1d), net(adj_idx1d));
                });
              }
          }
      }

      gen.allocate();

      for (std::size_t i{0} ; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; connected(l))
          if (pn_->inner_node(r)) { // macro-macro 
            auto [lr, rl] = gen.set(net(l), net(r));
            de_to_throat[lr] = i;
            de_to_throat[rl] = i;
            throat_to_de[i] = lr;
          }
          else if constexpr (Outlet) {
            if (pn_->outlet() == r) // macro-outlet
              gen.set(net(l), vertex_count - 1);
          }

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if constexpr (Outlet)
              if (voxel_idx_t adj_idx1d{map_idx(img_->dim.x() - 1, j, k)}; connected(adj_idx1d)) // darcy-outlet
                gen.set(net(adj_idx1d), vertex_count - 1);

            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (img_->velem[idx1d]) // macro-darcy 
                  gen.set(net(img_->velem[idx1d]), net(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx_t adj_idx1d = idx1d + map_idx[d]; img_->phase[adj_idx1d] == microporous) // darcy-darcy
                      gen.set(net(idx1d), net(adj_idx1d));
                });
              }
          }
      }

      return {gen.acquire(), std::move(de_to_throat), std::move(throat_to_de)};
    }

    std::tuple</*HYPRE_BigInt, */size_t, dpl::hypre::ls_known_storage> generate_pressure_input(
      idx1d_t nrows, const dpl::strong_vector<net_tag, idx1d_t>& forward, double const_permeability, auto filter) const {

      using namespace presets;


      auto voxel_idx_map = img_->idx1d_mapper();

      dpl::hypre::ls_known_storage_builder builder{
        nrows/**connected_count_*/,
        [&forward, this](auto i) { return forward[this->net(i)]; }
      };

      for (std::size_t i = 0; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i];
          pn_->inner_node(r) && connected(l) && filter(l) && filter(r) && filter(i)) // macro-macro // TODO: order
          builder.reserve(l, r);

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d) && filter(idx1d)) {
                if (auto velem = img_->velem[idx1d]; velem && filter(velem)) // macro-darcy
                  builder.reserve(velem, idx1d);

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (auto adj_idx1d = idx1d + voxel_idx_map[d]; img_->phase[adj_idx1d] == microporous && filter(adj_idx1d)) // darcy-darcy
                      builder.reserve(idx1d, adj_idx1d);
                });
              }
      }

      builder.allocate();

      for (std::size_t i = 0; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; connected(l) && filter(i) && filter(l)) {
          auto coef = pn_->coef(i);

          if (pn_->inner_node(r)) {// macro-macro
            if (filter(r))
              builder.set(l, r, coef);
          }
          else if (r == pn_->inlet()) { // macro-inlet      // TODO FILTER: INLET MUST BE FILLED
            builder.add_b(l, coef/**1 Pa*/);
            builder.add_diag(l, coef);
          }
          else { // macro-outlet                            // TODO FILTER: OUTLET MUST BE FILLED
            // builder.add_b(l, coef/**0 Pa*/);
            builder.add_diag(l, coef);
          }
        }

      {
        auto cell_size = pn_->physical_size/img_->dim;
        
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx_t idx1d{0};

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx_t inlet_idx1d = voxel_idx_map(0, j, k); connected(inlet_idx1d) && filter(inlet_idx1d)) { // darcy-inlet
              auto coef = -2*cell_size.x()*const_permeability;

              builder.add_b(inlet_idx1d, coef/**1 Pa*/);
              builder.add_diag(inlet_idx1d, coef);
            }

            if (voxel_idx_t outlet_idx1d = voxel_idx_map(img_->dim.x() - 1, j, k); connected(outlet_idx1d) && filter(outlet_idx1d)) { // darcy-outlet
              auto coef = -2*cell_size.x()*const_permeability;

              // builder.add_b(outlet_idx1d), coef/**0 Pa*/);
              builder.add_diag(outlet_idx1d, coef);
            }

            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d) && filter(idx1d)) {
                if (auto velem = img_->velem[idx1d]; velem && filter(velem)) { // macro-darcy
                  macro_idx_t adj_macro_idx = velem;

                  using eq_tri = hydraulic_properties::equilateral_triangle_properties;

                  auto li = cell_size.x()/2;
                  auto gi = const_permeability;
                  
                  auto lj = pn_->node_[attribs::r_ins][*adj_macro_idx];
                  auto gj = eq_tri::conductance(eq_tri::area(lj));
                  
                  auto lt = (cell_size*(ijk + 0.5) - pn_->node_[attribs::pos][*adj_macro_idx]).length() - li - lj;
                  auto gt = gj;
                  
                  if (lt < 0)
                    lt = 0;
                  
                  auto coef = -1.0/(li/gi + lt/gt + lj/gj);

                  builder.set(adj_macro_idx, idx1d, coef);
                }

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (auto adj_idx1d = idx1d + voxel_idx_map[d]; img_->phase[adj_idx1d] == microporous && filter(adj_idx1d)) { // darcy-darcy
                      auto coef = -cell_size.x()*const_permeability;
                      builder.set(idx1d, adj_idx1d, coef);
                    }
                });
              }
          }
      }

      return {/**connected_count_, */builder.nvalues(), builder.acquire()};
    }


   
    void flow_summary(const dpl::strong_vector<net_tag, HYPRE_Complex>& pressure, double const_permeability, auto filter) const {
      double inlet_flow_sum = 0;
      double outlet_flow_sum = 0;

      auto cell_size_x = (pn_->physical_size/img_->dim).x();

      auto map_idx = img_->idx1d_mapper();

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx_t inlet_idx1d = map_idx(0, j, k); connected(inlet_idx1d) && filter(inlet_idx1d)) // darcy-inlet
              inlet_flow_sum += -2*cell_size_x*const_permeability*(1 - pressure[net(inlet_idx1d)]);

            if (voxel_idx_t outlet_idx1d = map_idx(img_->dim.x() - 1, j, k); connected(outlet_idx1d) && filter(outlet_idx1d)) // darcy-outlet
              outlet_flow_sum += -2*cell_size_x*const_permeability*(pressure[net(outlet_idx1d)]);
          }
      }

      for (std::size_t i = 0; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; r == pn_->inlet()) { // macro-inlet
          if (connected(l) && filter(i) && filter(l))
            inlet_flow_sum += pn_->coef(i)*(1 - pressure[net(l)]);
        }
        else if (r == pn_->outlet()) { // macro-outlet
          if (connected(l) && filter(i) && filter(l))
            outlet_flow_sum += pn_->coef(i)*(pressure[net(l)]);
        }

      
      std::cout << fmt::format(
R"(microprs perm: {} mD
  inlet perm: {:.6f} mD
  outlet perm: {:.6f} mD
)",
        const_permeability/presets::darcy_to_m2*1000,
        -inlet_flow_sum/pn_->physical_size.x()/presets::darcy_to_m2*1000,
        -outlet_flow_sum/pn_->physical_size.x()/presets::darcy_to_m2*1000);
    }

    double eval_total_pore_volume(double darcy_porosity) const {
      auto volume = 0.0;

      for (macro_idx_t i{0}; i < pn_->node_count(); ++i)
        if (connected(i))
          volume += pn_->node_[attribs::volume][*i];

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && connected(l))
          volume += pn_->throat_[attribs::volume][i];
        
      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim).prod()*darcy_porosity;

      for (voxel_idx_t i{0}; i < img_->size; ++i)
        if (connected(i))
          volume += unit_darcy_pore_volume;

      return volume;
    }
  };
}
