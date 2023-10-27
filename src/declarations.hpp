#pragma once

#include <fstream>
#include <numbers>

#include <dpl/static_vector.hpp>
#include <dpl/hypre/mpi_module.hpp>

#include <HYPRE_utilities.h>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/pending/disjoint_sets.hpp>

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
  class pressure_cache
  {
    static inline std::unique_ptr<std::unordered_map<std::size_t, double>> cache_ = nullptr;

    static auto path() {
      return std::filesystem::path(dpl::hypre::mpi::mpi_exec).replace_filename("cache") / "solve.json";
    }

  public:
    static void load() {
      cache_ = std::make_unique<std::unordered_map<std::size_t, double>>();

      try {
        for (const auto& j : nlohmann::json::parse(std::ifstream{path()})) {
          (*cache_)[std::stoull(j["hash"].get<std::string>(), nullptr, 16)] = j["value"].get<double>();
        }
      }
      catch (...) {}
    }

    static void save() {
      nlohmann::json arr;

      for (auto [hash, value] : *cache_) {
        nlohmann::json j;
        j["hash"] = fmt::format("{:x}", hash);
        j["value"] = value;
        arr.push_back(j);
      }

      std::ofstream{path()} << arr.dump(2);
    }

    static auto& cache() {
      return *cache_;
    }
  };

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

    /**
     * \brief phase in the center
     */
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
    struct equilateral_triangle
    {
      static inline constexpr auto sqrt_3 = 1.732050807568877293527446341505872366942805253810380628055806; // std::sqrt(3)
      static inline constexpr auto beta = std::numbers::pi/3/2;
      static inline constexpr auto sin_beta = 0.5;
      static inline constexpr auto corners = 3;

      static constexpr double sqr(double x) { return x*x; }

      static constexpr double area(double r_ins = 1) { return 3*sqrt_3*sqr(r_ins); }
      static constexpr double perimeter(double r_ins = 1) { return 6*sqrt_3*r_ins; }
      static constexpr auto shape_factor() { return area()/sqr(perimeter()); }

      // k*G, k - coefficient, G - shape factor
      static double conductance_single(double area = 1, double viscosity = 1) {
        return sqrt_3/60*sqr(area)/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }

      static double conductance_films(double theta, double film_area = 1) {
        return sqr(film_area)*(
          0.364
          *(-0.261799387799 + 0.25*theta + 0.5*cos(theta)*cos(0.523598775598 + theta))
          /sqr(1.04719755120 - theta + 2.0*cos(0.523598775598 + theta))
          + 0.28*0.0481125224325
          )/3.0;
      }

      static constexpr bool has_films(double theta_rec) {
        return theta_rec < std::numbers::pi/3;
      }

      static double r_cap_piston_no_films(double theta, double r_ins = 1) {      
        return r_ins/(2*cos(theta));          
      }

      static double r_cap_snap_off(double theta, double r_ins = 1) {
        return r_ins/(cos(theta) - sin(theta)/sqrt_3);
      }

      static double area_films(double theta, double r_cap = 1) { // TODO: replace with area_corners_valv
        return sqr(r_cap)*(
          -0.5435164422364771 + 3*theta + 1.7320508075688772*std::cos(2*theta) +
          1.7320508075688772*std::sin(0.5235987755982988 - 2*theta)
        );
      }

      static double r_cap_piston_with_films(double theta, double r_ins = 1) {  // TODO: replace with Valvatne
        return r_ins/(std::cos(theta) + 0.759835685652*std::sqrt(1.04719755120 - theta + std::cos(theta)*std::sin(theta)));          
      }

      static double r_cap_piston_with_films_valvatne_full(double theta, double r_ins = 1) {
        using namespace std;

        auto cos_theta = cos(theta);

        auto S1 = corners*(cos_theta*cos(theta + beta)/sin_beta + theta + beta - numbers::pi/2);
        auto S2 = corners*(cos(theta + beta)/sin_beta);
        auto S3 = 2*corners*(numbers::pi/2 - theta - beta);

        auto D = S1 - 2*S2*cos_theta + S3;

        return r_ins/(1 + sqrt(1 + 4*shape_factor()*D/sqr(cos_theta)))/cos_theta;
        // return r_ins*cos_theta*(-1 + sqrt(1 + 4*G*D/sqr(cos_theta)))/(4*D*G);
      }

      static double r_cap_piston_with_films_valvatne(double theta, double r_ins = 1) {
        using namespace std;

        auto cos_theta = cos(theta);
        auto D = corners*(numbers::pi/2 - cos_theta*cos(theta + beta)/sin_beta - theta - beta);

        return r_ins/(1 + sqrt(1 + 4*shape_factor()*D/sqr(cos_theta)))/cos_theta;
      }

      static double b_length(double theta, double r_cap) { // ab length in promef
        return r_cap*cos(theta + beta)/sin_beta;
      }

      static double r_cap_hinging(double b, double theta_h) {
        return b*sin_beta/cos(theta_h + beta);
      }

      static double hinging(double b, double r_cap) {
        return std::acos(b*sin_beta/r_cap) - beta;
      }

      static double area_corners(double theta, double r_cap = 1) { // = area_films
        return sqr(r_cap)*corners*(cos(theta)*cos(theta + beta)/sin_beta + theta + beta - std::numbers::pi/2);
      }

      static double simple_balance(double theta, double r_cap, double r_ins) {
        using namespace std;
      
        auto A_eff = area(r_ins) - area_corners(theta, r_cap);
        auto L_os = r_ins/2/shape_factor() - 2*corners*b_length(theta, r_cap);
        auto L_ow = 2*r_cap*corners*(numbers::pi/2 - beta - theta);
      
        return A_eff/(L_ow + L_os*cos(theta)) - r_cap;
      }

      static double pinned_balance(double b_rec, double theta_adv, double r_cap, double r_ins) {
        using namespace std;

        auto theta_h = hinging(b_rec, r_cap);

        auto A_eff = area(r_ins) - area_corners(theta_h, r_cap);
        auto L_os = r_ins/2/shape_factor() - 2*corners*b_rec;
        auto L_ow = 2*r_cap*corners*asin(b_rec*sin_beta/r_cap);

        return A_eff/(L_ow + L_os*cos(theta_adv)) - r_cap;
      }

      /**
       * \brief
       *   theta_rec < Pi/3 (60 Deg)
       *   theta_rec <= theta_adv < Pi/3 (60 Deg)
       */
      static double r_cap_piston_secondary(double theta_rec, double r_cap_min, double theta_adv, double r_ins) {
        auto b_rec = b_length(theta_rec, r_cap_min);

        auto pc_for_hinge_adv = 1/r_cap_hinging(b_rec, theta_adv);
        auto pc_for_hinge_rec = 1/r_cap_min;

        auto func = [=](double pc) { return pinned_balance(b_rec, theta_adv, 1/pc, r_ins); };

        uintmax_t max = 64;

        auto fa = func(pc_for_hinge_adv);
        auto fb = func(pc_for_hinge_rec);

        if (fa*fb < 0) // hinging has solution
          return 1/boost::math::tools::toms748_solve(
            func,
            pc_for_hinge_adv,
            pc_for_hinge_rec,
            fa,
            fb,
            boost::math::tools::eps_tolerance<double>(),
            max).first;

        return r_cap_piston_with_films_valvatne(theta_adv, r_ins); // non-hinging
      }









      // static double r_cap_by_alpha(double ab, double alpha) {
      //   return ab/2.0/sin(alpha);
      // }
      //
      // static double alpha(double ab, double r_cap) { // loooks like alpha is not theta_h (hinging)
      //   /*
      //    * alpha = Pi/3 - thetaH
      //    * thetaH [0; 5Pi/6]
      //    * alpha [-Pi/2; Pi/3]
      //    *
      //    * |rCap| >= ab/2
      //   */
      //
      //   return asin(ab/2.0/r_cap);
      // }
      //
      // static double ab_length(double theta_rec, double r_cap) { // b_i in Valvatne
      //   return r_cap*(1.7320508075688772*cos(theta_rec) - sin(theta_rec));
      // }
      //
      // static double ab_length(double theta_rec, double ratio, double r_ins) {
      //   auto cos_theta = cos(theta_rec);
      //   auto sin_theta = sin(theta_rec);
      //   return r_ins*ratio*(1.7320508075688772*cos_theta - sin_theta)/
      //     (cos_theta + 0.7598356856515925*sqrt(1.0471975511965976 - theta_rec + cos_theta*sin_theta));
      // }
      //
      //
      // class pinned_to_single
      // {
      //   double ab;  
      //   double t0;
      //   double t1;
      //
      // public:
      //   // elem.drainageRatio = pnm.capillary_radius_global/equilateral_triangle_properties::rCap_PLD_with_wetting_films(elem.phase1ContactAngle[0], elem.rInsEffectivePc);
      //   // auto pinnedToSingle = phase_properties<DefPhase>::sign*equilateral_triangle_properties::rCap(elem.abLength, pinnedAlpha);
      //
      //   pinned_to_single(double theta_rec, double theta_adv, double ratio) {
      //     ab = ab_length(theta_rec, ratio);    
      //     t0 = (6.928203230275509- 4.*ab)*cos(theta_adv);
      //     t1 = 1.7320508075688772*(ab*ab - 4.0)/ab;
      //   }
      //
      //   std::tuple<double, double> operator()(double alpha) const {
      //     auto sinAlpha = sin(alpha);
      //     auto cosAlpha = cos(alpha);    
      //
      //     return std::make_tuple(
      //       t0 + t1*sinAlpha + ab*(cosAlpha + alpha/sinAlpha),
      //       t1*cosAlpha + ab*((1.0 - alpha*cosAlpha/sinAlpha)/sinAlpha - sinAlpha));
      //   }
      //
      //
      //   double find_alpha() const {
      //     using namespace std::numbers;
      //     return boost::math::tools::newton_raphson_iterate(*this, 1e-9, -2.0000*pi/3.0, 1.00000*pi/3.0, 22);  
      //   }
      // };
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






  






  // #define def_strong_type(name, value_type) \
  //   struct name##_tag {}; \
  //   using name##_t = dpl::strong_integer<value_type, name##_tag>;

  /**
   * \brief maximum node count
   */
  using idx1d_t = int32_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;

  struct voxel_tag {};
  using voxel_t = dpl::strong_integer<idx1d_t, voxel_tag>;

  struct macro_tag {};
  using macro_t = dpl::strong_integer<idx1d_t, macro_tag>;

  struct net_tag {};
  using net_t = dpl::strong_integer<idx1d_t, net_tag>;

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

      constexpr operator macro_t() const {
        return macro_t{value};
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

  
  namespace attrib {
    def_attrib(pos)
    def_attrib(r_ins)
    def_attrib(adj)
    def_attrib(length)
    def_attrib(length0)
    def_attrib(length1)
    def_attrib(volume)
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
    return map_idx3_t<voxel_t>{dim};
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
    bool save_cache = true;
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

    double darcy_perm = std::numeric_limits<double>::quiet_NaN(); /* mD */
    double darcy_poro = std::numeric_limits<double>::quiet_NaN(); /* fraction */

    struct input_curves
    {
      std::vector<dpl::vector2d> pc; /* [Sw, Pc] */
      std::array<std::vector<dpl::vector2d>, 2> kr;

      auto calc_pc_inv() const {
        using namespace std::ranges;
        auto inv = pc;
        reverse(inv);
        for_each(inv, [](dpl::vector2d& p) { std::swap(p.x(), p.y()); });
        inv.resize(unique(inv, {}, [](const dpl::vector2d& p) { return p.x(); }).begin() - inv.begin());
        return inv;
      }

    } primary,
      secondary;
    

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

      if (auto j_micro = j.find("microporosity"); j_micro != j.end()) {
        darcy_perm = (*j_micro)["permeability"].get<double>()*0.001*presets::darcy_to_m2;
        darcy_poro = (*j_micro)["porosity"];

        if (auto j_pr = j_micro->find("primary"); j_pr != j_micro->end()) {
          primary.pc = (*j_pr)["capillary_pressure"];
          primary.kr[0] = (*j_pr)["relative_permeability"][0];
          primary.kr[1] = (*j_pr)["relative_permeability"][1];  
        }

        if (auto j_sec = j_micro->find("secondary"); j_sec != j_micro->end()) {
          secondary.pc = (*j_sec)["capillary_pressure"];
          secondary.kr[0] = (*j_sec)["relative_permeability"][0];
          secondary.kr[1] = (*j_sec)["relative_permeability"][1];  
        }
      }

      solver.load(j["solver"]);
      loaded = true;
    }
  };

  inline void crop(
    const std::filesystem::path& src_p, const dpl::vector3i& src_size, const dpl::vector3i& src_origin,
    const std::filesystem::path& dst_p, const dpl::vector3i& dst_size) {

    auto src_total_size = src_size.prod();
    auto dst_total_size = dst_size.prod();

    std::vector<unsigned char> src(src_total_size);
    std::vector<unsigned char> dst(dst_total_size);

    std::ifstream is(src_p);
    
    is.read(reinterpret_cast<char*>(src.data()), src_total_size);

    auto src_mapper = idx_mapper(src_size);

    idx3d_t ijk;
    auto& [i, j, k] = ijk;
    voxel_t idx1d{0};

    for (k = 0; k < dst_size.z(); ++k)
      for (j = 0; j < dst_size.y(); ++j)
        for (i = 0; i < dst_size.x(); ++i, ++idx1d)
          dst[*idx1d] = src[*src_mapper(src_origin + ijk)];

    std::ofstream os(dst_p);
    os.write(reinterpret_cast<char*>(dst.data()), dst_total_size);
  }
}
