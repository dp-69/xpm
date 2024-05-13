#pragma once

#include <filesystem>
#include <fstream>
#include <numbers>
#include <regex>

#include <dpl/curve2d.hpp>
#include <dpl/fmt-formatter.hpp>
#include <dpl/json.hpp>
#include <dpl/hypre/core.hpp>
#include <dpl/qt/parse.hpp>

#include <boost/math/tools/roots.hpp>
#include <boost/pending/disjoint_sets.hpp>
  
  


namespace xpm
{
  class config_exception final : public std::runtime_error
  {
  public:
    // explicit config_exception(const std::string& _Message)
    //   : runtime_error(_Message) {}

    explicit config_exception(const char* msg)
      : runtime_error(msg) {}
  };

  class pressure_cache
  {
    static inline std::unique_ptr<std::unordered_map<std::size_t, double>> cache_ = nullptr;

    static auto path() {
      return std::filesystem::path(dpl::mpi::exec).replace_filename("cache")/"solve.json";
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

    static auto hash(std::size_t nvalues, const dpl::hypre::ls_known_storage& lks) {
      return std::hash<std::string_view>{}(
        std::string_view{
          reinterpret_cast<char*>(lks.values.get()),
          nvalues * sizeof(HYPRE_Complex)
        });
    }
  };

  class phase_config
  {
    static constexpr unsigned char phase0_ = 0;
    static constexpr unsigned char phase1_ = 128;

    static constexpr unsigned char phase_bits_ = 192;
    static constexpr unsigned char layout_bits_ = 3;
    

    static constexpr unsigned char single_ = 0;
    static constexpr unsigned char bulk_films_ = 1;
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




  

  namespace hydraulic_properties
  {
    struct equilateral_triangle
    {
      static constexpr auto sqrt_3 = 1.732050807568877293527446341505872366942805253810380628055806; // std::sqrt(3)
      static constexpr auto beta = std::numbers::pi/3/2;
      static constexpr auto sin_beta = 0.5;
      static constexpr auto corners = 3;

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
  using idx1d_t = int64_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;

  namespace helper
  {
    struct voxel_tag {};
    struct macro_tag {};
    struct net_tag {};
  }

  using voxel_t = dpl::strong_integer<idx1d_t, helper::voxel_tag>;
  using macro_t = dpl::strong_integer<idx1d_t, helper::macro_tag>;
  using net_t   = dpl::strong_integer<idx1d_t, helper::net_tag>;  

  using throat_t = std::size_t;

  template<class T>
  concept macro_voxel_t = std::is_same_v<T, macro_t> || std::is_same_v<T, voxel_t>;

  template<class T>
  concept macro_throat_t = std::is_same_v<T, macro_t> || std::is_same_v<T, throat_t>;

  
  

  
  

  namespace voxel_prop
  {
    namespace helper
    {
      struct phase_tag {};
      struct velem_tag {};
    }

    
    using phase_t = dpl::strong_integer<std::uint8_t, helper::phase_tag>;

    /**
     * \brief
     *   for a void voxel - a macro node it belongs
     *   for a microporous voxel - an adjacent pore node
     */
    struct velem_t : dpl::strong_integer<std::int32_t, helper::velem_tag, true>
    {
      constexpr velem_t() : strong_integer{invalid_value()} {}
      constexpr explicit velem_t(type v) : strong_integer{v} {}

      explicit constexpr operator macro_t() const {
        return macro_t{value};
      }
    };
  }

  namespace parse
  {
    struct image_dict
    {
      std::uint8_t pore;
      std::uint8_t solid;

      void load(const nlohmann::json& j) {
        pore = j["void"];
        solid = j["solid"];
      }

      bool is_void(const voxel_prop::phase_t p) const { return *p == pore; }
      bool is_solid(const voxel_prop::phase_t p) const { return *p == solid; }
      bool is_darcy(const voxel_prop::phase_t p) const { return !is_void(p) && !is_solid(p); }
    };
  }


  namespace presets {
    static constexpr auto darcy_to_m2 = 9.869233e-13;
    static constexpr auto mD_to_m2 = 9.869233e-16;
  }

  
  namespace attrib {
    def_attrib(pos)
    def_attrib(r_ins)

    /*
     * left index < right index
     *
     * inlet or outlet are always when (inlet, outlet) is not possible
     *
     */
    def_attrib(adj) 
    def_attrib(length)
    def_attrib(length0)
    def_attrib(length1)
    def_attrib(volume)
  }

  // using disjoint_sets = boost::disjoint_sets<int*, idx1d_t*>;
  //   boost::disjoint_sets_with_storage<
  //   boost::typed_identity_property_map<idx1d_t>,
  //   boost::typed_identity_property_map<idx1d_t>
  // >;

  template <typename Rank = idx1d_t>
  std::tuple<std::unique_ptr<Rank[]>, std::unique_ptr<idx1d_t[]>> init_rank_parent(idx1d_t size) {
    auto rank = std::make_unique<Rank[]>(size);
    auto parent = std::make_unique<idx1d_t[]>(size);

    std::fill_n(rank.get(), size, 0);
    std::iota(parent.get(), parent.get() + size, 0);

    return std::make_tuple(std::move(rank), std::move(parent));
  }

  template <typename>
  struct wrapper {};

  template <>
  struct wrapper<nlohmann::json>
  {
    const nlohmann::json* json = nullptr;

    wrapper() = default;
    wrapper(const nlohmann::json* j) : json(j) {}
    wrapper(const nlohmann::json& j) : json(&j) {}

    wrapper operator()(const auto& key, const auto&... rest) const {
      if (auto jj = (*this)(key); jj)
        return jj(rest...);
      return {};
    }

    auto operator()(const auto& key) const {
      auto found = json->find(key);
      return found == json->end() ? wrapper{} : wrapper{&*found};
    }

    auto& operator*() const {
      return *json;
    }

    auto operator->() const {
      return json;
    }

    explicit operator bool() const {
      return json;
    }

    // void operator>>(auto& value) {
    //   if (json)
    //     value = *json;
    // }

    // wrapper find(std::string_view path) {
    //   auto pos = path.find('/');
    //   while (pos != std::string_view::npos) {
    //     path.remove_prefix(pos + 1);
    //     return (*this)(std::string_view(path.data(), pos)).find(path);
    //   }
    //   return (*this)(path);
    // }

    auto set(auto& value, const auto&... keys) {
      if (auto jj = (*this)(keys...); jj)
        value = *jj;
    }
  };

  // void operator<<(auto& value, const wrapper<nlohmann::json>& arg) {
  //   if (arg)
  //     value = *arg;
  // }

  // void operator<<(auto& value, const wrapper<nlohmann::json>& arg) {
  //   if (arg)
  //     value = *arg;
  // }

  struct poro_perm_t
  {
    double poro;
    double perm;
    QColor color;


    static poro_perm_t nan() {
      poro_perm_t pp;
      pp.poro = std::numeric_limits<double>::quiet_NaN();
      pp.perm = std::numeric_limits<double>::quiet_NaN();
      return pp;
      // return poro_perm_t{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    bool is_nan() const {
      return std::isnan(poro);
    }
  };

  struct runtime_settings
  {
    using wrap = wrapper<nlohmann::json>;

    bool loaded = false;
    bool occupancy_images = false;

    double max_pc = std::numeric_limits<double>::max();

    struct {
      std::filesystem::path path;
      dpl::vector3i size;
      double resolution;
      parse::image_dict phases;
      // bool grey = false;

      void load(const nlohmann::json& j) {
        path = std::string{j["path"]};
        size = j["size"];
        resolution = j["resolution"];
        phases.load(j["phase"]);
        // wrap{j}.set(grey, "grey");
      }

      auto pnextract_filename() {
        auto stem = path.stem();

        if (std::regex_match(stem.string(), std::regex{R"(.+_(\d+)x(\d+)x(\d+)_(\d+)p(\d+)um$)"}))
          return path.filename();

        auto resol = fmt::format("{:.3f}", resolution/1e-6);
        resol[resol.find('.')] = 'p';

        return std::filesystem::path{
          fmt::format("{}_{:d}x{:d}x{:d}_{}um{}",
            stem, size.x(), size.y(), size.z(), resol, path.extension())
        };
      }
    } image;


    double theta = 0;

    struct
    {
      // double perm_single = std::numeric_limits<double>::quiet_NaN();

      // double A = std::numeric_limits<double>::quiet_NaN();
      // double n1 = 1;
      // double n2 = 1;

      dpl::strong_array<voxel_prop::phase_t, poro_perm_t> poro_perm{poro_perm_t::nan()};

      void set_poro_perm(const nlohmann::json& list) {
        using namespace voxel_prop;

        int i = 0;

        auto mult = 255./(list.size() + 1);  // NOLINT(clang-diagnostic-implicit-int-float-conversion)

        for (const auto& j : list) {
          auto& ref = poro_perm[phase_t{j["value"].get<phase_t::type>()}];

          ref.poro = j["porosity"];
          ref.perm = j["permeability"].get<double>()*presets::mD_to_m2;
          ref.color = QColor::fromHsl(++i*mult, 175, 122);  // NOLINT(cppcoreguidelines-narrowing-conversions)

          dpl::qt::try_parse(j, ref.color);

          // poro_perm[phase_t{j["value"].get<phase_t::type>()}] =
          //   {j["porosity"], j["permeability"].get<double>()*presets::mD_to_m2};
        }
      }

      auto poro(voxel_prop::phase_t p) const {
        return poro_perm[p].poro;
      }

      auto perm(voxel_prop::phase_t p) const {
        return poro_perm[p].perm;
      }
    } darcy;

    
    struct input_curves {
      dpl::curve2d pc;  /* [Sw, Pc] */
      std::array<dpl::curve2d, 2> kr;
    } primary,
      secondary;
    

    struct {
      std::optional<dpl::vector3i> decomposition;
      HYPRE_Real tolerance = 1.e-20;
      HYPRE_Int max_iterations = 20;
      HYPRE_Int aggressive_levels = 0;

      struct {
        bool use = true;
        bool save = true;
      } cache;

      void load(wrap j) {
        tolerance = (*j)["tolerance"];
        max_iterations = (*j)["max_iterations"];
        j.set(aggressive_levels, "aggressive_number_of_levels");
        j.set(decomposition, "decomposition");
        j.set(cache.use, "cache", "use");
        j.set(cache.save, "cache", "save");
      }
    } solver;


    struct {
      bool invasion_percolation = true;
      std::string display = "saturation";
      double sw_of_pc = 0.05;
      double sw_of_kr = 0.075;
    } report;
    

    void load(wrap j) {
      image.load(*j("image"));

      if (auto j_micro = j("microporosity"); j_micro) {
        // const auto& first_record = (*j_micro)["voxel"][0];
        // darcy.perm_single = first_record["permeability"].get<double>()*0.001*presets::darcy_to_m2;
        // darcy.poro_single = first_record["porosity"];

        darcy.set_poro_perm((*j_micro)["voxel"]);

        // j_micro.set(darcy.n1, "kozeny_carman", "n1");
        // j_micro.set(darcy.n2, "kozeny_carman", "n2");
        // darcy.A = darcy.perm_single*std::pow(1 - darcy.poro_single, darcy.n2)/std::pow(darcy.poro_single, darcy.n1);


        if (auto j_primary = j_micro("primary"); j_primary) {
          primary.pc = {(*j_primary)["capillary_pressure"]};
          primary.kr[0] = {(*j_primary)["relative_permeability"][0]};
          primary.kr[1] = {(*j_primary)["relative_permeability"][1]};  
        }

        if (auto j_secondary = j_micro("secondary"); j_secondary) {
          secondary.pc = {(*j_secondary)["capillary_pressure"]};
          secondary.kr[0] = {(*j_secondary)["relative_permeability"][0]};
          secondary.kr[1] = {(*j_secondary)["relative_permeability"][1]};  
        }
      }

      j.set(occupancy_images, "report", "occupancy_images");
      j.set(report.sw_of_pc, "report", "capillary_pressure_sw_step");
      j.set(report.sw_of_kr, "report", "relative_permeability_sw_step");
      j.set(report.display, "report", "display");
      j.set(report.invasion_percolation, "report", "invasion_percolation");
      j.set(max_pc, "max_capillary_pressure");

      if (auto j_theta = j("macro_contact_angle"); j_theta)
        theta = j_theta->get<double>()/180*std::numbers::pi;

      solver.load(*j("solver"));
      loaded = true;
    }
  };

  template<typename Type = std::uint8_t, typename Proj = std::identity>
  void transform(
    const std::filesystem::path& src_path, const dpl::vector3i& src_size, const dpl::vector3i& src_origin,
    const std::filesystem::path& dst_path, const dpl::vector3i& dst_size, Proj proj = {}) {

    auto src_total_size = src_size.prod();
    auto dst_total_size = dst_size.prod();

    std::vector<Type> src(src_total_size);
    std::vector<Type> dst(dst_total_size);

    
    
    std::ifstream is{src_path, std::ios::binary};
    
    is.read(reinterpret_cast<char*>(src.data()), src_total_size*sizeof(Type));

    dpl::idx1d_map<voxel_t> src_mapper{src_size};

    idx3d_t ijk;
    auto& [i, j, k] = ijk;
    voxel_t idx1d{0};

    for (k = 0; k < dst_size.z(); ++k)
      for (j = 0; j < dst_size.y(); ++j)
        for (i = 0; i < dst_size.x(); ++i, ++idx1d)
          dst[*idx1d] = proj(src[*src_mapper(src_origin + ijk)]);

    std::ofstream{dst_path, std::ios::binary}
      .write(reinterpret_cast<char*>(dst.data()), dst_total_size*sizeof(Type));
  }
}
