#pragma once

#include <dpl/static_vector.hpp>

#include <boost/iostreams/device/mapped_file.hpp>

#include <iostream>
#include <fstream>

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

      const auto& operator*() const {
        return value;
      }
    };

    /**
     * \brief
     *   -1  : for non-pore voxels
     *   >=0 : pore node of a pore-voxel
     */
    struct velem
    {
      std::int32_t value;

      // auto& operator*() {
      //   return value;
      // }

      const auto& operator*() const {
        return value;
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


  /**
   * \brief maximum node count
   */
  using idx1d_t = int32_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;



  //TODO
  struct idx1d_expl
  {
    idx1d_t value;
  };

  //TODO
  struct idx3d_expl
  {
    idx3d_t value;
  };

  // template<typename T, int n> requires std::integral<T>
  // class map_idx_t
  // {
  //   // dpl::vector_n<Type, n> dim;
  // };

  template <typename R> requires std::integral<R>
  class map_idx3_t
  {
    R x_, xy_;

  public:
    template<typename T>
    explicit map_idx3_t(dpl::vector_n<T, 3> dim)
      : x_(dim.x()), xy_(static_cast<R>(dim.x())*dim.y()) {}

    template<typename V>
    R operator()(dpl::vector_n<V, 3> v) const {
      return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    }

    template<typename V>
    R operator()(V x, V y, V z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    template <typename I>
    auto operator[](std::integral_constant<I, 0>) { return 1; }
     
    template<typename I>
    auto operator[](std::integral_constant<I, 1>) { return x_; }

    template <typename I>
    auto operator[](std::integral_constant<I, 2>) { return xy_; }
  };


  inline auto idx_mapper(idx3d_t dim) {
    return map_idx3_t<idx1d_t>{dim};
  }














  struct image_data
  {
    idx3d_t dim;
    idx1d_t size;

    std::unique_ptr<voxel_tag::phase[]> phase;
    /**
     * \brief
     *    for a pore voxel, velem is a pore node it belongs
     *    for a microporous voxel, velem is an adjacent pore node
     */
    std::unique_ptr<voxel_tag::velem[]> velem;

    auto idx1d_mapper() const {
      return idx_mapper(dim);
    }

    void eval_microporous_velem() {
      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      idx1d_t idx1d = 0;

      auto map_idx = idx1d_mapper();

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i, ++idx1d) {
            int32_t adj = -1;
                
            if (phase[idx1d] == presets::microporous) {
              dpl::sfor<3>([&](auto d) {
                if (ijk[d] > 0)
                  if (auto adj_idx = idx1d - map_idx[d]; phase[adj_idx] == presets::pore)
                    if (auto adj_velem = *velem[adj_idx]; adj_velem >= 0)
                      adj = adj_velem;

                if (ijk[d] < dim[d] - 1)
                  if (auto adj_idx = idx1d + map_idx[d]; phase[adj_idx] == presets::pore)
                    if (auto adj_velem = *velem[adj_idx]; adj_velem >= 0)
                      adj = adj_velem;
              });

              velem[idx1d] = {adj};
            }
          }

      #ifdef XPM_DEBUG_OUTPUT  
        std::cout << "\n\nVelems_adj array produced and filled"; // NOLINT(clang-diagnostic-misleading-indentation)
      #endif
    }
  };


  


  namespace parse
  {
    struct image_dict
    {
      std::uint8_t solid;
      std::uint8_t pore;
      std::uint8_t microporous;
    };

    std::pair<std::unique_ptr<voxel_tag::phase[]>, std::size_t> read_image(const auto& image_path, image_dict input_config) {
      std::ifstream is(image_path);
      is.seekg(0, std::ios_base::end);
      size_t image_size = is.tellg();
      is.seekg(0, std::ios_base::beg);
      auto phase_arr = std::make_unique<voxel_tag::phase[]>(image_size);
      is.read(reinterpret_cast<char*>(phase_arr.get()), image_size);

      size_t pore_voxels = 0;
      size_t solid_voxels = 0;
      size_t microporous_voxels = 0;

      for (size_t i = 0; i < image_size; ++i)
        if (auto& value = phase_arr[i];
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

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << std::format("\n\ntotal: {}; pore: {}; solid: {}; microporous: {} voxels", image_size, pore_voxels, solid_voxels, microporous_voxels);
      #endif

      return {std::move(phase_arr), image_size};
    }

    /**
     * \brief
     * input file value description
     *   -2: solid (validated),
     *   -1: inlet/outlet (do not know?),
     *   0, 1: do not exist (validated),
     *   >=2: cluster (0, 1 are inlet and outlet node indices in Statoil format)
     */
    inline std::unique_ptr<voxel_tag::velem[]> read_icl_velems(const std::filesystem::path& network_path, const idx3d_t& dim) {
      boost::iostreams::mapped_file_source file(network_path.string() + "_VElems.raw");
      const auto* file_ptr = reinterpret_cast<const std::int32_t*>(file.data());

      auto result = std::make_unique<voxel_tag::velem[]>(dim.prod());
      auto* ptr = result.get();

      idx3d_t velems_factor{1, dim.x() + 2, (dim.x() + 2)*(dim.y() + 2)};
      idx3d_t ijk;
      auto& [i, j, k] = ijk;

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j) 
          for (i = 0; i < dim.x(); ++i) {
            auto val = file_ptr[velems_factor.dot(ijk + 1)];
            *ptr++ = {val > 0 ? val - 2 : -1/*val*/};
          }

      return result;
    }
  }
}