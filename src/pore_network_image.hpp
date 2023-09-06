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
  struct row_decomposition
  {
    std::unique_ptr<idx1d_t[]> net_to_decomposed;  // TODO: strong_array
    std::unique_ptr<idx1d_t[]> decomposed_to_net;
    std::vector<dpl::hypre::index_range> blocks;
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
    
    
    idx1d_t parse_statoil_text_idx(std::integral auto idx) const {
      if (idx == -1)
        return node_count();
  
      if (idx == 0)
        return node_count() + 1;
    
      return idx - 1;
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
      attribs::adj_t, std::pair<macro_idx, macro_idx>,     
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

    


    
    
    idx1d_t node_count() const {
      return static_cast<idx1d_t>(node_.size());
    }

    size_t throat_count() const {
      return throat_.size();
    }

    bool inner_node(macro_idx i) const {
      return i < node_count();
    }

    auto inlet() const {
      return macro_idx{node_count()};
    }

    auto outlet() const {
      return macro_idx{node_count() + 1};
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
          *left = parse_statoil_text_idx(value); 

          parse_text(throat1_ptr, value);
          *right = parse_statoil_text_idx(value);        
          
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


    double coef(size_t i) const {
      using eq_tri = hydraulic_properties::equilateral_triangle_properties;
      using namespace attribs;

      auto [left, right] = throat_[adj][i];

      return -1.0/(
        throat_[length0][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][*left])) +
        throat_[length][i]/eq_tri::conductance(eq_tri::area(throat_[r_ins][i])) +
        (inner_node(right) ? throat_[length1][i]/eq_tri::conductance(eq_tri::area(node_[r_ins][*right])) : 0.0));
    }

    bool eval_inlet_outlet_connectivity() const {
      disjoint_sets ds(node_count() + 2);

      for (auto [l, r] : throat_.range(attribs::adj))
        ds.union_set(*l, *r);

      return ds.find_set(*inlet()) == ds.find_set(*outlet());
    }

    std::pair<dpl::hypre::ls_known_storage, size_t> generate_pressure_input() const {
      dpl::hypre::ls_known_storage_builder builder;

      builder.allocate_rows(node_count());

      for (auto i : dpl::range(throat_count()))
        if (auto [l, r] = throat_[attribs::adj][i]; inner_node(r))
          builder.reserve_connection(*l, *r);

      builder.allocate_values();

      for (auto i : dpl::range(throat_count())) {
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
          builder.set_connection(*l, *r, coef);
      }

      

      return {builder.acquire_storage(), builder.nvalues()};
    }


     
    void connectivity_flow_summary(const dpl::hypre::mpi::solve_result& solve  /*const std::unique_ptr<double[]>& pressure*/) const {
      auto& [pressure, residual, iter] = solve; 

      disjoint_sets ds(node_count());
      std::vector<bool> connected_inlet(node_count());
      std::vector<bool> connected_outlet(node_count());
      
      auto connected = [&](macro_idx i) {
        auto rep = ds.find_set(*i);
        return connected_inlet[rep] && connected_outlet[rep];
      };
      
      {
        for (auto [l, r] : throat_.range(attribs::adj))
          if (inner_node(r))
            ds.union_set(*l, *r);
      
        for (auto [l, r] : throat_.range(attribs::adj))
          if (r == inlet())
            connected_inlet[ds.find_set(*l)] = true;
          else if (r == outlet())
            connected_outlet[ds.find_set(*l)] = true;
      
        idx1d_t disconnected_macro = 0;
        
        for (macro_idx i{0}; i < node_count(); ++i)
          if (!connected(i))
            ++disconnected_macro;
      
        double inlet_flow = 0;
        double outlet_flow = 0;
      
        std::cout << fmt::format("  isolated {} nodes", disconnected_macro);
      
        for (auto i : dpl::range(throat_count())) {
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
      auto nrows = node_count();

      dpl::hypre::mpi::save(input, nrows, nvalues, {{0, nrows - 1}}, tolerance, max_iterations);

      std::system(fmt::format("mpiexec -np 1 \"{}\" -s",  // NOLINT(concurrency-mt-unsafe)
        dpl::hypre::mpi::mpi_exec).c_str());

      connectivity_flow_summary(dpl::hypre::mpi::load_values(nrows));
    }


    void connectivity_flow_summary(HYPRE_Real tolerace, HYPRE_Int max_iterations) const {
      auto pressure = std::make_unique<double[]>(node_count());
      auto [residual, iters] = dpl::hypre::solve(
        {0, node_count() - 1}, generate_pressure_input().first, pressure.get(), tolerace, max_iterations); // gross solve (with isolated)
      connectivity_flow_summary({std::move(pressure), residual, iters});
    }
  };




  struct image_data
  {
    idx3d_t dim;
    idx1d_t size;

    std::unique_ptr<voxel_property::phase[]> phase;

    /**
     * \brief
     *    for a void voxel - a macro node it belongs
     *    for a microporous voxel - an adjacent pore node
     */
    std::unique_ptr<voxel_property::velem[]> velem;

    auto idx1d_mapper() const {
      return idx_mapper(dim);
    }

    void eval_microporous() {
      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      voxel_idx idx1d{0};

      auto map_idx = idx1d_mapper();

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i, ++idx1d) {
            int32_t adj = -1;
                
            if (phase[*idx1d] == presets::microporous) {
              dpl::sfor<3>([&](auto d) {
                if (ijk[d] > 0)
                  if (voxel_idx adj_idx = idx1d - map_idx[d]; phase[*adj_idx] == presets::pore)
                    if (auto adj_velem = *velem[*adj_idx]; adj_velem >= 0)
                      adj = adj_velem;

                if (ijk[d] < dim[d] - 1)
                  if (voxel_idx adj_idx = idx1d + map_idx[d]; phase[*adj_idx] == presets::pore)
                    if (auto adj_velem = *velem[*adj_idx]; adj_velem >= 0)
                      adj = adj_velem;
              });

              *velem[*idx1d] = adj;
            }
          }
    }

    void read_image(const auto& image_path, parse::image_dict input_config) {
      std::ifstream is(image_path);
      is.seekg(0, std::ios_base::end);
      size = static_cast<idx1d_t>(is.tellg());
      is.seekg(0, std::ios_base::beg);
      phase = std::make_unique<voxel_property::phase[]>(size);
      is.read(reinterpret_cast<char*>(phase.get()), size);

      size_t pore_voxels = 0;
      size_t solid_voxels = 0;
      size_t microporous_voxels = 0;

      for (auto i : dpl::range(size))
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

      velem = std::make_unique<voxel_property::velem[]>(dim.prod());

      auto* ptr = velem.get();

      idx3d_t velems_factor{1, dim.x() + 2, (dim.x() + 2)*(dim.y() + 2)};
      idx3d_t ijk;
      auto& [i, j, k] = ijk;

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i) {
            auto val = file_ptr[velems_factor.dot(ijk + 1)];
            **ptr++ = val > 0 ? val - 2 : -1/*val*/;
          }
    }
  };



  class pore_network_image
  {
    pore_network* pn_;
    image_data* img_;

    static inline constexpr net_idx isolated_idx_{std::numeric_limits<idx1d_t>::max()};

    /**
     * \brief
     *   .size() = pn_->node_count_ + img_->size
     *   maps total() index to compressed effective (flowing) index
     */
    dpl::strong_array<total_tag, net_tag, idx1d_t> net_map_;

    /**
     * \brief effective porosity elements count
     */
    net_idx connected_count_{0};


    

  public:
    void init(pore_network* pn, image_data* img) {
      pn_ = pn;
      img_ = img;
    }

    // idx1d_t inlet() const { return pn_->node_count() + img_->size_microporous; }
    // idx1d_t outlet() const { return inlet() + 1; }

    auto connected_count() const { return connected_count_; }

    auto total(macro_idx i) const { return total_idx{*i}; }
    auto total(voxel_idx i) const { return total_idx{pn_->node_count() + *i}; }

    auto net(macro_idx i) const { return net_map_[total(i)]; }
    auto net(voxel_idx i) const { return net_map_[total(i)]; }

    bool connected(macro_idx i) const {
      return net(i) != isolated_idx_;
    }

    /**
     * \brief if connected, the voxel must be microporous, i.e. the latter evaluation is redundant
     */
    bool connected(voxel_idx i) const {
      return net(i) != isolated_idx_;
    }

    void evaluate_isolated() {
      auto gross_total_size = pn_->node_count() + img_->size;

      disjoint_sets ds(gross_total_size);

      for (auto [l, r] : pn_->throat_.range(attribs::adj))
        if (pn_->inner_node(r)) // macro-macro
          ds.union_set(*l, *r);

      using namespace presets;
      auto map_idx = img_->idx1d_mapper();

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (img_->phase[*idx1d] == microporous) {
                if (macro_idx adj_macro_idx = img_->velem[*idx1d]; adj_macro_idx >= 0) // macro-darcy
                  ds.union_set(*total(adj_macro_idx), *total(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx adj_idx1d = idx1d + map_idx[d]; img_->phase[*adj_idx1d] == microporous) // darcy-darcy
                      ds.union_set(*total(idx1d), *total(adj_idx1d));
                });
              }
      }

      std::vector<bool> inlet(gross_total_size);
      std::vector<bool> outlet(gross_total_size);

      for (auto [l, r] : pn_->throat_.range(attribs::adj))
        if (r == pn_->inlet()) // macro-inlet
          inlet[ds.find_set(*total(l))] = true;
        else if (r == pn_->outlet()) // macro-outlet
          outlet[ds.find_set(*total(l))] = true;

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx inlet_idx1d{map_idx(0, j, k)}; img_->phase[*inlet_idx1d] == microporous) // darcy-inlet
              inlet[ds.find_set(*total(inlet_idx1d))] = true;

            if (voxel_idx outlet_idx1d{map_idx(img_->dim.x() - 1, j, k)}; img_->phase[*outlet_idx1d] == microporous) // darcy-outlet
              outlet[ds.find_set(*total(outlet_idx1d))] = true;
          }
      }

      net_map_.resize(gross_total_size);

      for (auto i : dpl::range(gross_total_size)) {
        auto rep = ds.find_set(i); 
        net_map_[total_idx{i}] = inlet[rep] && outlet[rep] ? connected_count_++ : isolated_idx_;
      }
    }


    row_decomposition decompose_rows(v3i blocks) const {
      auto block_size = pn_->physical_size/blocks;

      using pair = std::pair<idx1d_t, int>;

      std::vector<pair> net_idx_block(*connected_count_);

      auto map_idx = idx_mapper(blocks);

      idx1d_t net_idx = 0;

      auto add = [&](const v3d& pos) {
        v3i b_idx = pos/block_size;

        net_idx_block[net_idx] = {net_idx,
          *map_idx(
            std::clamp(b_idx.x(), 0, blocks.x() - 1), 
            std::clamp(b_idx.y(), 0, blocks.y() - 1), 
            std::clamp(b_idx.z(), 0, blocks.z() - 1))
        };

        ++net_idx;
      };

      {
        for (macro_idx i{0}; i < pn_->node_count(); ++i)
          if (connected(i)) // macro node
            add(pn_->node_[attribs::pos][*i]);


        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        using namespace presets;

        auto cell_size = pn_->physical_size/img_->dim;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) // darcy node
                add(cell_size*(ijk + 0.5));
      }

      std::ranges::sort(net_idx_block, [](const pair& l, const pair& r) { return l.second < r.second; });

      row_decomposition mapping;

      mapping.net_to_decomposed = std::make_unique<idx1d_t[]>(*connected_count_);
      mapping.decomposed_to_net = std::make_unique<idx1d_t[]>(*connected_count_);

      auto block_count = blocks.prod();
      auto row_count_per_block = std::make_unique<idx1d_t[]>(block_count);
      mapping.blocks.resize(block_count);

      for (auto i : dpl::range(block_count))
        row_count_per_block[i] = 0;
      
      for (auto i : dpl::range(*connected_count_)) {
        mapping.net_to_decomposed[net_idx_block[i].first] = i;
        mapping.decomposed_to_net[i] = net_idx_block[i].first;
        ++row_count_per_block[net_idx_block[i].second];
      }

      HYPRE_BigInt first_row = 0;
      for (auto i : dpl::range(block_count)) {
        mapping.blocks[i] = {first_row, first_row + row_count_per_block[i] - 1};
        first_row += row_count_per_block[i];
      }

      return mapping;
    }


    dpl::graph::dc_graph generate_dc_graph() const {
      using namespace dpl::graph;
      using namespace presets;
      auto map_idx = img_->idx1d_mapper();

      net_idx vertex_count = connected_count_ + 1; // last is outlet
      net_idx outlet_idx = vertex_count - 1;

      graph_generator<net_idx> gen(*vertex_count);

      for (auto [l, r] : pn_->throat_.range(attribs::adj)/* | std::ranges::views::filter([](auto pair) { return true; })*/)
        if (connected(l))
          if (pn_->inner_node(r)) // macro-macro
            gen.reserve(net(l), net(r));
          else if (pn_->outlet() == r) // macro-outlet
            gen.reserve(net(l), outlet_idx);


      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx adj_idx1d{map_idx(img_->dim.x() - 1, j, k)}; connected(adj_idx1d)) // darcy-outlet
              gen.reserve(net(adj_idx1d), outlet_idx);

            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (macro_idx adj_macro_idx{img_->velem[*idx1d]}; adj_macro_idx >= 0) // macro-darcy 
                  gen.reserve(net(adj_macro_idx), net(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx adj_idx1d = idx1d + map_idx[d]; img_->phase[*adj_idx1d] == microporous) // darcy-darcy
                      gen.reserve(net(idx1d), net(adj_idx1d));
                });
              }
          }
      }

      gen.allocate();

      for (auto [l, r] : pn_->throat_.range(attribs::adj))
        if (connected(l))
          if (pn_->inner_node(r)) // macro-macro
            gen.set(net(l), net(r));
          else if (pn_->outlet() == r) // macro-outlet
            gen.set(net(l), outlet_idx);

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx adj_idx1d{map_idx(img_->dim.x() - 1, j, k)}; connected(adj_idx1d)) // darcy-outlet
              gen.set(net(adj_idx1d), outlet_idx);

            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (macro_idx adj_macro_idx = img_->velem[*idx1d]; adj_macro_idx >= 0) // macro-darcy 
                  gen.set(net(adj_macro_idx), net(idx1d));

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx adj_idx1d = idx1d + map_idx[d]; img_->phase[*adj_idx1d] == microporous) // darcy-darcy
                      gen.set(net(idx1d), net(adj_idx1d));
                });
              }
          }
      }

      return gen.acquire();
    }

    std::tuple<HYPRE_BigInt, size_t, dpl::hypre::ls_known_storage> generate_pressure_input(
      const row_decomposition& mapping, double const_permeability) const {

      using namespace attribs;
      using namespace presets;

      auto map_idx = img_->idx1d_mapper();

      auto* block = mapping.net_to_decomposed.get();
      
      dpl::hypre::ls_known_storage_builder builder;

      builder.allocate_rows(*connected_count_);
      
      for (auto [l, r] : pn_->throat_.range(adj))
        if (pn_->inner_node(r) && connected(l)) // macro-macro
          builder.reserve_connection(
            block[*net(l)],
            block[*net(r)]);

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j)
            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (macro_idx adj_macro_idx = img_->velem[*idx1d]; adj_macro_idx >= 0) // macro-darcy
                  builder.reserve_connection(
                    block[*net(adj_macro_idx)],
                    block[*net(idx1d)]);

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx adj_idx1d = idx1d + map_idx[d]; img_->phase[*adj_idx1d] == microporous) // darcy-darcy
                      builder.reserve_connection(
                        block[*net(idx1d)],
                        block[*net(adj_idx1d)]);
                });
              }
      }


      builder.allocate_values();


      for (auto i : dpl::range(pn_->throat_count()))
        if (auto [l, r] = pn_->throat_[adj][i]; connected(l)) {
          auto coef = pn_->coef(i);

          if (pn_->inner_node(r)) // macro-macro
            builder.set_connection(
              block[*net(l)],
              block[*net(r)],
              coef);
          else if (r == pn_->inlet()) { // macro-inlet
            builder.add_b(block[*net(l)], coef/**1 Pa*/);
            builder.add_diag(block[*net(l)], coef);
          }
          else { // macro-outlet
            // builder.add_b(block[*net(l)], coef/**0 Pa*/);
            builder.add_diag(block[*net(l)], coef);
          }
        }


      {
        auto cell_size = pn_->physical_size/img_->dim;
        
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_idx idx1d;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx inlet_idx1d = map_idx(0, j, k); connected(inlet_idx1d)) { // darcy-inlet
              auto coef = -2*cell_size.x()*const_permeability;

              builder.add_b(block[*net(inlet_idx1d)], coef/**1 Pa*/);
              builder.add_diag(block[*net(inlet_idx1d)], coef);
            }

            if (voxel_idx outlet_idx1d = map_idx(img_->dim.x() - 1, j, k); connected(outlet_idx1d)) { // darcy-outlet
              auto coef = -2*cell_size.x()*const_permeability;

              // builder.add_b(block[*net(outlet_idx1d)], coef/**0 Pa*/);
              builder.add_diag(block[*net(outlet_idx1d)], coef);
            }


            for (i = 0; i < img_->dim.x(); ++i, ++idx1d)
              if (connected(idx1d)) {
                if (macro_idx adj_macro_idx = img_->velem[*idx1d]; adj_macro_idx >= 0) { // macro-darcy
                  using eq_tri = hydraulic_properties::equilateral_triangle_properties;

                  // ReSharper disable CppInconsistentNaming
                  auto Li = cell_size.x()/2;
                  auto gi = const_permeability;
                  
                  auto Lj = pn_->node_[r_ins][*adj_macro_idx];
                  auto gj = eq_tri::conductance(eq_tri::area(Lj));
                  
                  auto Lt = (cell_size*(ijk + 0.5) - pn_->node_[pos][*adj_macro_idx]).length() - Li - Lj;
                  auto gt = gj;
                  // ReSharper restore CppInconsistentNaming
                  
                  if (Lt < 0)
                    Lt = 0;
                    // std::cout << "NEGATIVE Lt < 0";
                  
                  auto coef = -1.0/(Li/gi + Lt/gt + Lj/gj);

                  // auto length = (cell_size*(ijk + 0.5) - pn_->node_[pos][*adj_macro_idx]).length(); // NOLINT(clang-diagnostic-shadow)
                  // auto r_ins = cell_size.x()/4;                                                    // NOLINT(clang-diagnostic-shadow)
                  // auto coef = -1.0/(length/eq_tri::conductance(eq_tri::area(r_ins)));

                  builder.set_connection(
                    block[*net(adj_macro_idx)],
                    block[*net(idx1d)],
                    coef);
                }

                dpl::sfor<3>([&](auto d) {
                  if (ijk[d] < img_->dim[d] - 1)
                    if (voxel_idx adj_idx1d = idx1d + map_idx[d]; img_->phase[*adj_idx1d] == microporous) { // darcy-darcy
                      auto coef = -cell_size.x()*const_permeability;

                      builder.set_connection(
                        block[*net(idx1d)],
                        block[*net(adj_idx1d)],
                        coef);
                  }
                });
              }
          }
      }

      return {*connected_count_, builder.nvalues(), builder.acquire_storage()};
    }


   
    void flow_summary(const double* pressure, double const_permeability) const {
      double inlet_flow_sum = 0;
      double outlet_flow_sum = 0;

      auto cell_size_x = (pn_->physical_size/img_->dim).x();

      using namespace presets;
      auto map_idx = img_->idx1d_mapper();

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;

        for (k = 0; k < img_->dim.z(); ++k)
          for (j = 0; j < img_->dim.y(); ++j) {
            if (voxel_idx inlet_idx1d = map_idx(0, j, k); connected(inlet_idx1d)) // darcy-inlet
              inlet_flow_sum += -2*cell_size_x*const_permeability*(1 - pressure[*net(inlet_idx1d)]);

            if (voxel_idx outlet_idx1d = map_idx(img_->dim.x() - 1, j, k); connected(outlet_idx1d)) // darcy-outlet
              outlet_flow_sum += -2*cell_size_x*const_permeability*(pressure[*net(outlet_idx1d)]);
          }
      }

      for (auto i : dpl::range(pn_->throat_count()))
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; r == pn_->inlet()) { // macro-inlet
          if (connected(l))
            inlet_flow_sum += pn_->coef(i)*(1 - pressure[*net(l)]);
        }
        else if (r == pn_->outlet()) { // macro-outlet
          if (connected(l))
            outlet_flow_sum += pn_->coef(i)*(pressure[*net(l)]);
        }

      
      std::cout << fmt::format(
R"(microprs perm: {} mD
  inlet perm: {:.6f} mD
  outlet perm: {:.6f} mD
)",
        const_permeability/darcy_to_m2*1000,
        -inlet_flow_sum/pn_->physical_size.x()/darcy_to_m2*1000,
        -outlet_flow_sum/pn_->physical_size.x()/darcy_to_m2*1000);
    }
  };
}
