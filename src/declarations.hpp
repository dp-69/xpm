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

template <typename T, typename Tag>
struct fmt::formatter<dpl::strong_integer<T, Tag>> : formatter<T>
{
  template <typename FormatContext>
  auto format(dpl::strong_integer<T, Tag> si, FormatContext& ctx) {
    return formatter<T>::format(*si, ctx);
  }
};

  
  


namespace xpm
{

  class phase_config
  {
    static inline constexpr unsigned char phase0_ = 0;
    static inline constexpr unsigned char phase1_ = 128;

    static inline constexpr unsigned char phase_bits_ = 192;
    static inline constexpr unsigned char layout_bits_ = 3;
    

    static inline constexpr unsigned char single_ = 0;
    static inline constexpr unsigned char bulk_films_ = 1;
    // static inline constexpr unsigned char phase1_bulk_films_ = phase1_ | bulk_films_;

    unsigned char value;

  public:

    constexpr explicit phase_config(unsigned char v = 0) : value(v) {}

    auto phase() const {
      return value & phase_bits_;
    }

    auto layout() const {
      return value & layout_bits_;
    }

    static constexpr auto phase0() {
      return phase0_;
    }

    static constexpr auto phase1() {
      return phase1_;
    }

    static constexpr auto bulk_films() {
      return bulk_films_;
    }

    static constexpr auto phase1_bulk_films() {
      return phase_config{phase1_ | bulk_films_};
    }
  };


  /**
   * \brief
   *   (00)00000(00)
   *    |       | layout
   *    |
   *    | bulk phase
   */
  // enum class phase_config : unsigned char
  // {
  //   phase0 = 0, 
  //   phase1 = 128,
  //
  //   single = 0,
  //   bulk_films = 1, // phase0 default
  //
  //   phase1_bulk_films = phase1 | bulk_films
  //
  //   // bulk = 1,
  //   // layers = 2, // only defending phase, intermediate stage
  //   // corners = 4,
  //   //
  //   // bulk_corners = bulk | corners, // 5
  //   // bulk_layers_corners = bulk | layers | corners, // 7
  //   //
  //   // b = bulk,
  //   // bc = bulk_corners,
  //   // blc = bulk_layers_corners
  // };

  // constexpr phase_config operator|(phase_config l, phase_config r) {
  //   using t = std::underlying_type_t<phase_config>;
  //   return static_cast<phase_config>(static_cast<t>(l) | static_cast<t>(r));
  // }




  class default_maps
  {
    static inline auto true_ = [](auto) { return true; };
    static inline auto unity_ = [](auto) { return 1; };

  public:
    using true_t = decltype(true_);
    using unity_t = decltype(unity_);

    static bool invert(std::true_type, bool v) { return !v; }
    static bool invert(std::false_type, bool v) { return v; }
  };

  namespace hydraulic_properties
  {
    class equilateral_triangle_properties
    {
      static double sqr(double x) { return x*x; }

    public:
      // struct shape
      // {
      //   static constexpr double area(double r_ins = 1) {
      //     return 5.19615242271*r_ins*r_ins;
      //   }
      // };

      static double area(double r_ins = 1) {
        return 5.19615242271*sqr(r_ins);
      }

      // k*G, k - coefficient, G - shape factor
      static double conductance_single_phase(double area = 1, double viscosity = 1) {
        return 0.0288675134595*sqr(area)/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }

      static double conductance_films(double theta, double film_area = 1) {
        return sqr(film_area)*(
          0.364
          *(-0.261799387799 + 0.25*theta + 0.5*cos(theta)*cos(0.523598775598 + theta))
          /sqr(1.04719755120 - theta + 2.0*cos(0.523598775598 + theta))
          + 0.28*0.0481125224325
          )/3.0;
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
        return sqr(r_cap)*(
          -0.5435164422364771 + 3*theta + 1.7320508075688772*std::cos(2*theta) +
          1.7320508075688772*std::sin(0.5235987755982988 - 2*theta)
        );
      }
    };
  }


  // props.conductanceSinglePhaseCoef = equilateral_triangle_properties::conductance_single_phase_coef()/**sqr(props.total_area())*/;
  //       
  //     props.conductanceWettingCornerCoef = 
  //        (0.364*(-0.261799387799 + 0.25*theta + 0.5*cos(theta)*cos(0.523598775598 + theta))/
  //              sqr(1.04719755120 - theta + 2.0*cos(0.523598775598 + theta)) +
  //        0.28*0.0481125224325)/3.0;   

  // props.wettingAreaSaturation_rCapSqr = equilateral_triangle_properties::total_corner_area_uniform_wet(theta)/props.total_area();

  //return props.phases.bulk_phase == Phase
  //        ? props.conductanceSinglePhaseCoef*(1 - props.wettingAreaSaturation_rCapSqr*sqr(rCap))          
  //        : props.conductanceWettingCornerCoef*sqr(props.wettingAreaSaturation_rCapSqr*sqr(rCap)); 






  









  /**
   * \brief maximum node count
   */
  using idx1d_t = int32_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;

  struct voxel_tag {};
  using voxel_idx_t = dpl::strong_integer<idx1d_t, voxel_tag>;

  struct macro_tag {};
  using macro_idx_t = dpl::strong_integer<idx1d_t, macro_tag>;

  struct net_tag {};
  using net_idx_t = dpl::strong_integer<idx1d_t, net_tag>;

  namespace voxel_property
  {
    struct phase_tag {};
    using phase_t = dpl::strong_integer<std::uint8_t, phase_tag>;

    /**
     * \brief
     *   -1  : for a non-void voxel
     *   >=0 : macro node of a void voxel
     */
    struct velem_tag {};
    struct velem_t : dpl::strong_integer<std::int32_t, velem_tag>
    {
    private:
      static inline constexpr auto not_valid = std::numeric_limits<value_type>::max();

    public:
      velem_t() : strong_integer(not_valid) {}
      constexpr explicit velem_t(const value_type v) : strong_integer(v) {}

      constexpr operator macro_idx_t() const {
        return macro_idx_t{value};
      }

      explicit constexpr operator bool() const {
        return value != not_valid;
      }
    };
  }

  namespace presets
  {
    static inline constexpr auto pore = voxel_property::phase_t{0};
    static inline constexpr auto solid = voxel_property::phase_t{1};
    static inline constexpr auto microporous = voxel_property::phase_t{2};

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

  template <typename R>
  class map_idx3_t
  {
    R x_, xy_;

  public:
    map_idx3_t() = default;

    template <typename T>
    explicit map_idx3_t(dpl::vector_n<T, 3> dim)
      : x_(dim.x()), xy_(static_cast<R>(dim.x())*dim.y()) {}

    // template <typename V>
    // R operator()(dpl::vector_n<V, 3> v) const {
    //   return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    // }

    template <typename T>
    R operator()(T x, T y, T z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    template <typename T>
    R operator()(const dpl::vector_n<T, 3>& v) const {
      return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    }

    auto operator()(std::integral_constant<int, 0>) const { return 1; }
    auto operator()(std::integral_constant<int, 1>) const { return x_; }
    auto operator()(std::integral_constant<int, 2>) const { return xy_; }
  };

  inline auto idx_mapper(idx3d_t dim) {
    return map_idx3_t<voxel_idx_t>{dim};
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
    bool use_cache = true;
    bool save_cache = false;
    bool loaded = false;

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

    double micro_perm = -999; /* mD */
    double micro_poro = -999; /* fraction */
    std::vector<dpl::vector2d> micro_pc;
    std::vector<dpl::vector2d> micro_kr0;
    std::vector<dpl::vector2d> micro_kr1;
    

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

      if (auto microprs = j.find("microporosity"); microprs != j.end()) {
        micro_perm = (*microprs)["permeability"];
        micro_poro = (*microprs)["porosity"];
        micro_pc = (*microprs)["capillary_pressure"];
        micro_kr0 = (*microprs)["relative_permeability"][0];
        micro_kr1 = (*microprs)["relative_permeability"][1];
      }

      solver.load(j["solver"]);
      loaded = true;
    }
  };

  inline void crop(
    const std::filesystem::path& src_p, dpl::vector3i src_size, dpl::vector3i src_origin,
    const std::filesystem::path& dst_p, dpl::vector3i dst_size) {

    auto src_total_size = src_size.prod();
    auto dst_total_size = dst_size.prod();

    std::vector<unsigned char> src(src_total_size);
    std::vector<unsigned char> dst(dst_total_size);

    std::ifstream is(src_p);
    
    is.read(reinterpret_cast<char*>(src.data()), src_total_size);

    auto src_mapper = idx_mapper(src_size);

    idx3d_t ijk;
    auto& [i, j, k] = ijk;
    voxel_idx_t idx1d{0};

    for (k = 0; k < dst_size.z(); ++k)
      for (j = 0; j < dst_size.y(); ++j)
        for (i = 0; i < dst_size.x(); ++i, ++idx1d)
          dst[*idx1d] = src[*src_mapper(src_origin + ijk)];

    std::ofstream os(dst_p);
    os.write(reinterpret_cast<char*>(dst.data()), dst_total_size);


    // if (connected(idx1d) && filter(idx1d)) {
    //   if (auto velem = img_->velem[idx1d]; velem && filter(velem)) // macro-darcy
    //     builder.reserve(velem, idx1d);
    //
    //   dpl::sfor<3>([&](auto d) {
    //     if (ijk[d] < img_->dim()[d] - 1)
    //       if (auto adj_idx1d = idx1d + img_->idx_map(d); img_->phase[adj_idx1d] == presets::microporous && filter(adj_idx1d)) // darcy-darcy
    //         builder.reserve(idx1d, adj_idx1d);
    //   });
    // }
  }
}
