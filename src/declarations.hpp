#pragma once

#include <dpl/static_vector.hpp>

#include <HYPRE_utilities.h>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <numbers>
#include <iostream>
#include <fstream>

#include <fmt/format.h>

template <>
struct fmt::formatter<std::filesystem::path> : formatter<std::string_view>
{
  template <typename FormatContext>
  auto format(const std::filesystem::path& path, FormatContext& ctx) {
    return formatter<std::string_view>::format(path.string(), ctx);
  }
};

namespace xpm
{
  namespace hydraulic_properties
  {
    struct equilateral_triangle_properties
    {
      // struct shape
      // {
      //   static constexpr double area(double r_ins = 1) {
      //     return 5.19615242271*r_ins*r_ins;
      //   }
      // };

      static double area(double r_ins = 1) {
        return 5.19615242271*r_ins*r_ins;
      }

      // k*G, k - coefficient, G - shape factor
      static double conductance(double area = 1, double viscosity = 1) {
        return 0.0288675134595*area*area/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }

      /**
       * \brief checks existence of wetting films
       */
      static bool has_films(double theta) {
        return theta < std::numbers::pi/3;
      }

      static double r_cap_piston_with_films(double theta, double r_ins = 1) {      
        return r_ins/(std::cos(theta) + 0.759835685652*std::sqrt(1.04719755120 - theta + std::cos(theta)*std::sin(theta)));          
      }

      static double area_of_films(double theta, double r_cap = 1) {
        return r_cap*r_cap*(
          -0.5435164422364771 + 3.*theta + 1.7320508075688772*std::cos(2.*theta) +
          1.7320508075688772*std::sin(0.5235987755982988 - 2.*theta)
        );
      }
    };
  }








  









  using v3i = dpl::vector3i;
  using v3d = dpl::vector3d;

  /**
   * \brief maximum node count
   */
  using idx1d_t = int32_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;

  struct voxel_tag {};
  using voxel_idx = dpl::strong_integer<idx1d_t, voxel_tag>;

  struct macro_tag {};
  using macro_idx = dpl::strong_integer<idx1d_t, macro_tag>;

  struct total_tag {};
  using total_idx = dpl::strong_integer<idx1d_t, total_tag>;

  struct net_tag {};
  using net_idx = dpl::strong_integer<idx1d_t, net_tag>;

  namespace voxel_property
  {
    struct phase_t {};
    using phase = dpl::strong_integer<std::uint8_t, phase_t>;

    /**
     * \brief
     *   -1  : for a non-void voxel
     *   >=0 : macro node of a void voxel
     */
    struct velem_t {};
    struct velem : dpl::strong_integer<std::int32_t, velem_t>
    {
      constexpr operator macro_idx() const { return macro_idx{value}; }
    };
  }

  namespace presets
  {
    static inline constexpr auto pore = voxel_property::phase{0};
    static inline constexpr auto solid = voxel_property::phase{1};
    static inline constexpr auto microporous = voxel_property::phase{2};

    static constexpr auto darcy_to_m2 = 9.869233e-13;
  }

  namespace attribs {
    def_static_key(pos)
    def_static_key(r_ins)
    def_static_key(adj)
    def_static_key(length)
    def_static_key(length0)
    def_static_key(length1)
    def_static_key(volume)
  }

  using disjoint_sets = boost::disjoint_sets_with_storage<
    boost::typed_identity_property_map<idx1d_t>,
    boost::typed_identity_property_map<idx1d_t>
  >;

  template <typename R>// requires std::integral<R>
  class map_idx3_t
  {
    R x_, xy_;

  public:
    template <typename T>
    explicit map_idx3_t(dpl::vector_n<T, 3> dim)
      : x_(dim.x()), xy_(static_cast<R>(dim.x())*dim.y()) {}

    template <typename V>
    R operator()(dpl::vector_n<V, 3> v) const {
      return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    }

    template <typename V>
    R operator()(V x, V y, V z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    auto operator[](std::integral_constant<int, 0>) { return 1; }
    auto operator[](std::integral_constant<int, 1>) { return x_; }
    auto operator[](std::integral_constant<int, 2>) { return xy_; }
  };

  inline auto idx_mapper(idx3d_t dim) {
    return map_idx3_t<voxel_idx>{dim};
  }
  
  namespace parse
  {
    struct image_dict
    {
      std::uint8_t pore;
      std::uint8_t solid;
      std::uint8_t microporous;

      void load(const nlohmann::json& j) {
        pore = j["void"];
        solid = j["solid"];
        microporous = j["microporous"];
      }
    };
  }

  struct startup_settings
  {
    struct {
      std::filesystem::path path;
      dpl::vector3i size;
      double resolution;
      parse::image_dict phases;

      void load(const nlohmann::json& j) {
        path = static_cast<std::string>(j["path"]);
        size = j["size"];
        resolution = j["resolution"];
        phases.load(j["phase"]);
      }
    } image;

    double microporous_perm; /* mD */

    struct
    {
      std::optional<dpl::vector3i> decomposition;
      HYPRE_Real tolerance = 1.e-20;
      HYPRE_Int max_iterations = 20;

      void load(const nlohmann::json& j) {
        tolerance = j["tolerance"];
        max_iterations = j["max_iterations"];

        if (auto d = j.find("decomposition"); d != j.end())
          decomposition = *d;
      }
    } solver;

    void load(const nlohmann::json& j) {
      image.load(j["image"]);
      microporous_perm = j["properties"]["microporosity"]["permeability"].get<double>();
      solver.load(j["solver"]);
    }
  };
}
