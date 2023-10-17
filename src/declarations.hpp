#pragma once

#include <dpl/static_vector.hpp>

#include <HYPRE_utilities.h>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <numbers>
#include <iostream>
#include <fstream>
#include <qmath.h>

#include <fmt/format.h>

#include <boost/math/tools/roots.hpp>

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
    struct equilateral_triangle_properties
    {
      static constexpr auto sqrt_3 = 1.732050807568877293527446341505872366942805253810380628055806; // std::sqrt(3)
      static constexpr auto beta = std::numbers::pi/3/2;
      static constexpr auto sin_beta = 0.5;
      static constexpr auto corners = 3;

      static constexpr double sqr(double x) { return x*x; }

      static constexpr double area(double r_ins = 1) {
        return 3*sqrt_3*sqr(r_ins);
      }

      static constexpr double perimeter(double r_ins = 1) {
        return 6*sqrt_3*r_ins;
      }

      static constexpr auto shape_factor() {
        return area()/sqr(perimeter());
      }

      // k*G, k - coefficient, G - shape factor
      static double conductance_single_phase(double area = 1, double viscosity = 1) {
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

      // static double r_cap_collapse(double theta_adv, double r_ins = 1) {
      //   return r_ins*1.73205080757/(1.73205080757*cos(theta_adv) - sin(theta_adv));
      // }

      // static double r_cap_piston_no_films_valvatne(double theta, double r_ins = 1) {
      //   static constexpr auto sqrt_pi_G = 0.3887800753890535401856827809473435907176788021745640816548206175;
      //   return r_ins/std::cos(theta)/(1 + 2*sqrt_pi_G);
      // }

      static constexpr bool has_films(double theta_rec) {
        return theta_rec < std::numbers::pi/3;
      }

      static double r_cap_piston_no_films(double theta, double r_ins = 1) {      
        return r_ins/(2*cos(theta));          
      }

      static double r_cap_snap_off_valv(double theta, double r_ins = 1) {
        return r_ins/(cos(theta) - sin(theta)/sqrt_3);
      }

      static double area_of_films(double theta, double r_cap = 1) { // TODO: replace with area_corners_valv
        return sqr(r_cap)*(
          -0.5435164422364771 + 3*theta + 1.7320508075688772*std::cos(2*theta) +
          1.7320508075688772*std::sin(0.5235987755982988 - 2*theta)
        );
      }

      static double r_cap_piston_with_films(double theta_rec, double r_ins = 1) {  // TODO: replace with Valvatne
        return r_ins/(std::cos(theta_rec) + 0.759835685652*std::sqrt(1.04719755120 - theta_rec + std::cos(theta_rec)*std::sin(theta_rec)));          
      }

      static double r_cap_piston_with_films_valvatne(double theta, double r_ins = 1) {
        using namespace std;

        auto cos_theta = cos(theta);

        auto S1 = corners*(cos_theta*cos(theta + beta)/sin_beta + theta + beta - numbers::pi/2);
        auto S2 = corners*(cos(theta + beta)/sin_beta);
        auto S3 = 2*corners*(numbers::pi/2 - theta - beta);

        auto D = S1 - 2*S2*cos_theta + S3;

        // return r_ins*cos_theta*(-1 + sqrt(1 + 4*G*D/sqr(cos_theta)))/(4*D*G);

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

      static double area_corners_valv(double theta, double r_cap = 1) { // = area_of_films
        return sqr(r_cap)*corners*(cos(theta)*cos(theta + beta)/sin_beta + theta + beta - std::numbers::pi/2);
      }

      static double simple_balance(double theta, double r_cap, double r_ins) {
        using namespace std;
      
        auto A_eff = area(r_ins) - area_corners_valv(theta, r_cap);
        auto L_os = r_ins/2/shape_factor() - 2*corners*b_length(theta, r_cap);
        auto L_ow = 2*r_cap*corners*(numbers::pi/2 - beta - theta);
      
        return A_eff/(L_ow + L_os*cos(theta)) - r_cap;
      }

      static double pinned_balance(double b_rec, double theta_adv, double r_cap, double r_ins) {
        using namespace std;

        auto theta_h = hinging(b_rec, r_cap);

        auto A_eff = area(r_ins) - area_corners_valv(theta_h, r_cap);
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

    double darcy_perm = -999; /* mD */
    double darcy_poro = -999; /* fraction */
    std::vector<dpl::vector2d> darcy_pc; /* [Sw, Pc] */
    std::vector<dpl::vector2d> darcy_kr0;
    std::vector<dpl::vector2d> darcy_kr1;
    

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

      if (auto micro = j.find("microporosity"); micro != j.end()) {
        darcy_perm = (*micro)["permeability"].get<double>()*0.001*presets::darcy_to_m2;
        darcy_poro = (*micro)["porosity"];
        darcy_pc = (*micro)["capillary_pressure"];
        darcy_kr0 = (*micro)["relative_permeability"][0];
        darcy_kr1 = (*micro)["relative_permeability"][1];
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
