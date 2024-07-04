/*
 * This file is part of Extensive Pore Modelling (xpm).
 *   | https://github.com/dp-69/xpm
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * xpm is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * xpm is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with xpm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <dpl/curve2d.hpp>
#include <dpl/defs_fmt.hpp>
#include <dpl/defs_json.hpp>
#include <dpl/hypre/core.hpp>
#include <dpl/qt/parse.hpp>

#include <boost/math/tools/roots.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <filesystem>
#include <fstream>
#include <numbers>
#include <regex>
  

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
          nvalues*sizeof(HYPRE_Complex)
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
      static constexpr double sqrt_3 = 1.732050807568877293527446341505872366942805253810380628055806; // std::sqrt(3)
      static constexpr double beta = std::numbers::pi/3/2;
      static constexpr double sin_beta = 0.5;
      static constexpr    int corners = 3;

      static constexpr double sqr(double x) { return x*x; }

      static constexpr double area(double r_ins = 1) { return 3*sqrt_3*sqr(r_ins); }
      static constexpr double perimeter(double r_ins = 1) { return 6*sqrt_3*r_ins; }
      static constexpr double shape_factor() { return area()/sqr(perimeter()); }

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

      // static double r_cap_piston_with_films(double theta, double r_ins = 1) {  // TODO: replace with Valvatne
      //   return r_ins/(std::cos(theta) + 0.759835685652*std::sqrt(1.04719755120 - theta + std::cos(theta)*std::sin(theta)));          
      // }

      // static double r_cap_piston_with_films_valvatne_full(double theta, double r_ins = 1) {
      //   using namespace std;
      //
      //   auto cos_theta = cos(theta);
      //
      //   auto S1 = corners*(cos_theta*cos(theta + beta)/sin_beta + theta + beta - numbers::pi/2);
      //   auto S2 = corners*(cos(theta + beta)/sin_beta);
      //   auto S3 = 2*corners*(numbers::pi/2 - theta - beta);
      //
      //   auto D = S1 - 2*S2*cos_theta + S3;
      //
      //   return r_ins/(1 + sqrt(1 + 4*shape_factor()*D/sqr(cos_theta)))/cos_theta;
      //   // return r_ins*cos_theta*(-1 + sqrt(1 + 4*G*D/sqr(cos_theta)))/(4*D*G);
      // }

      static double r_cap_piston_with_films_valvatne(double theta, double r_ins = 1) {
        auto cos_theta = cos(theta);
        auto D = corners*(std::numbers::pi/2 - cos_theta*cos(theta + beta)/sin_beta - theta - beta);

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
    struct voxel_tag { using type = idx1d_t; };
    struct macro_tag { using type = idx1d_t; };
    struct net_tag   { using type = idx1d_t; };
  }

  using voxel_t = dpl::so_integer<helper::voxel_tag>;
  using macro_t = dpl::so_integer<helper::macro_tag>;
  using net_t   = dpl::so_integer<helper::net_tag>;  

  using throat_t = std::size_t;

  template<typename T>
  concept macro_voxel_t = std::is_same_v<T, macro_t> || std::is_same_v<T, voxel_t>;

  template<typename T>
  concept macro_throat_t = std::is_same_v<T, macro_t> || std::is_same_v<T, throat_t>;

  
  
  
  
  

  namespace voxel_ns
  {
    namespace helper
    {
      struct phase_tag {
        using type = std::uint8_t;
      };

      struct velem_tag {
        using type = std::int32_t;
        static constexpr auto invalid_value = std::numeric_limits<type>::max();
      };
    }
    
    using phase_t = dpl::so_integer<helper::phase_tag>;
    static inline constexpr dpl::full_range_unsigned<phase_t> phases;

    /**
     * \brief
     *   for a void voxel  - a macro node it belongs to
     *   for a darcy voxel - an adjacent macro node or else <invalid> 
     */
    using velem_t = dpl::so_integer<helper::velem_tag>;

    // struct velem_t : dpl::so_integer<helper::velem_tag>
    // {
    //   velem_t() = default;
    //   constexpr explicit velem_t(const type v) : so_integer(v) {}
    //
    //   // explicit constexpr operator macro_t() const {
    //   //   return macro_t{value};
    //   // }
    //
    //   static constexpr auto invalid() {
    //     return velem_t{helper::velem_tag::invalid_value /*invalid_value()*/};
    //   }
    // };
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

  template <typename T> requires (std::is_same_v<std::remove_const_t<T>, nlohmann::json>)
  struct wrapper<T> 
  {
    T* json = nullptr;

    wrapper() = default;
    wrapper(T* j) : json(j) {}   // NOLINT(CppNonExplicitConvertingConstructor)
    wrapper(T& j) : json(&j) {}  // NOLINT(CppNonExplicitConvertingConstructor)

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

    auto try_set(auto& value, const auto&... keys) {
      if (auto jj = (*this)(keys...); jj)
        value = *jj;
    }
  };

  struct darcy_info
  {
    double poro;
    double perm;
    QColor color;

    std::array<dpl::curve2d, 2> pc_to_sw;
    std::array<std::array<dpl::curve2d, 2>, 2> kr;

    darcy_info() = default;

    darcy_info(const darcy_info& other) = delete;
    darcy_info(darcy_info&& other) noexcept = delete;
    darcy_info& operator=(const darcy_info& other) = delete;
    darcy_info& operator=(darcy_info&& other) noexcept = delete;
  };

  class runtime_settings
  {
    using json = nlohmann::json;

    static json load(const std::filesystem::path& filename, bool ignore_comments) {
      return json::parse(std::ifstream{filename}, nullptr, true, ignore_comments);
    }

    static const json& try_var(json& vars, const json& j, bool ignore_comments = false) {
      if (j.is_string()) {
        if (auto iter = vars.find(j); iter != vars.end())
          return *iter;

        return vars[j] = load(j, ignore_comments);
      }

      if (j.is_array() && j.size() == 2 && j[0].is_string() && j[1].is_string())
        return vars[j.dump()] = load(j[0], ignore_comments).at(j[1]);

      return j;
    }

    static json get_vars(wrapper<json> j, bool ignore_comments = false) {
      auto vars = j("variables");

      if (!vars)
        return json{};

      for (auto& v : *vars)
        if (v.is_string())
          v = load(v, ignore_comments);
        else if (v.is_array() && v.size() == 2 && v[0].is_string() && v[1].is_string())
          v = load(v[0], ignore_comments).at(v[1]);

      return std::move(*vars);
    }

  public:
    bool occupancy_images = false;

    double max_pc = std::numeric_limits<double>::max();
    double macro_mult = 1.0;

    struct image_settings {
      using phase_t = voxel_ns::phase_t;

      std::filesystem::path path;
      dpl::vector3i size;
      double resolution;

      void load(const json& j) {
        path = std::string{j["path"]};
        size = j["size"];
        resolution = j["resolution"];
        void_v  = j["phase"]["void"];
        solid_v = j["phase"]["solid"];
      }

      std::string pnextract_modelname() const {
        auto stem_str = path.stem().string();

        if (std::smatch matches;
          std::regex_match(stem_str, matches, std::regex{R"((.+)_(\d+)x(\d+)x(\d+)_(\d+)p(\d+)um$)"}))
          return matches[1];

        return stem_str;
      }

      auto pnextract_filename() const {
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

      /*
       * wide
       */
      phase_t void_v;
      phase_t solid_v;

      struct {                                              // NOLINT(cppcoreguidelines-pro-type-member-init)
        /*
         * wide
         */
        dpl::so_array<phase_t, bool> found;    
        dpl::so_array<phase_t, phase_t> narrow;

        /*
         * narrow
         */
        phase_t count;
        dpl::so_uptr<phase_t, darcy_info> info;

        auto span() const {
          return info.span(count);
        }
      } darcy;


      void set_poro_perm(json& list, json& vars) {
        auto s = list.size();

        if (s > ((1 << sizeof(phase_t)*8) - 2))
          throw config_exception("reached maximum number of darcy values.");

        darcy.count = phase_t(s);                          // NOLINT(clang-diagnostic-implicit-int-conversion)
        darcy.info.resize(darcy.count);    

        auto mult = 255./(s + 1);                          // NOLINT(cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-int-float-conversion)

        using namespace dpl;
        using namespace std;

        unordered_map<string, json> cache;

        for (phase_t i{0}; auto& j : list) {
          auto end = j.end();

          if (auto f = j.find("file"); f != end) {
            json* file_json;  // NOLINT(cppcoreguidelines-init-variables)

            string filename = *f;
            if (auto iter = cache.find(filename); iter == cache.end()) {
              file_json = &cache[filename];
              *file_json = json::parse(ifstream{filename});
            }
            else
              file_json = &iter->second;

            for (const auto& [key, val] : file_json->items())
              if (j.find(key) == end)
                j[key] = val;
          }

          auto p = j.at("value").get<phase_t>();

          darcy.narrow[p] = i;
          darcy.found[p] = true;

          auto& info = darcy.info[i];  // NOLINT(CppUseStructuredBinding)

          info.poro = j.at("poro");
          info.perm = j.at("perm").get<double>()*presets::mD_to_m2;

          if (auto iter = j.find("cap_press"); iter != end)
            for (int c = 0; c < 2; ++c){
              parse(try_var(vars, (*iter)[c]), info.pc_to_sw[c]);
              info.pc_to_sw[c] = info.pc_to_sw[c].inverse();
            }

          if (auto iter = j.find("rel_perm"); iter != end)
            for (int c = 0; c < 2; ++c) {
              const auto& rel_perm = try_var(vars, (*iter)[c]);
              parse(rel_perm.at(0), info.kr[c][0]);
              parse(rel_perm.at(1), info.kr[c][1]);
            }

          {
            info.color = QColor::fromHsl(*i*mult, 175, 122);  // NOLINT(cppcoreguidelines-narrowing-conversions)
            qt::try_parse(j, info.color);
          }

          ++i;
        }

        darcy.narrow[void_v] = darcy.count;
        *darcy.narrow[solid_v] = *darcy.count + 1;
      }
    } image;

    double theta = 0;

    struct {
      std::optional<dpl::vector3i> decomposition;
      HYPRE_Real tolerance = 1.e-20;
      HYPRE_Int max_iterations = 20;
      HYPRE_Int aggressive_levels = 0;
      HYPRE_Int print_level = 0;

      struct {
        bool use = true;
        bool save = true;
      } cache;

      void load(wrapper<const json> j) {
        tolerance = (*j)["tolerance"];
        max_iterations = (*j)["max_iterations"];
        j.try_set(aggressive_levels, "aggressive_number_of_levels");
        j.try_set(print_level, "print_level");
        j.try_set(decomposition, "decomposition");
        j.try_set(cache.use, "cache", "use");
        j.try_set(cache.save, "cache", "save");
      }
    } solver;


    struct {
      bool invasion_percolation = true;
      std::string display = "saturation";
      double sw_of_pc = 0.05;
      double sw_of_kr = 0.075;
    } report;


    void load(wrapper<json> j) {
      image.load(*j("image"));

      auto vars = get_vars(j);

      if (auto jj = j("darcy"); jj) {
        image.set_poro_perm(*jj, vars);

        // j_micro.set(darcy.n1, "kozeny_carman", "n1");
        // j_micro.set(darcy.n2, "kozeny_carman", "n2");
        // darcy.A = darcy.perm_single*std::pow(1 - darcy.poro_single, darcy.n2)/std::pow(darcy.poro_single, darcy.n1);
      }
      else {
        image.darcy.count = voxel_ns::phase_t{0};
        image.darcy.info.resize(voxel_ns::phase_t{0});
      }

      j.try_set(occupancy_images, "report", "occupancy_images");
      j.try_set(report.sw_of_pc, "report", "capillary_pressure_sw_step");
      j.try_set(report.sw_of_kr, "report", "relative_permeability_sw_step");
      j.try_set(report.display, "report", "display");
      j.try_set(report.invasion_percolation, "report", "invasion_percolation");
      j.try_set(max_pc, "max_capillary_pressure");
      j.try_set(macro_mult, "macro", "trans_multiplier");

      if (auto j_theta = j("macro_contact_angle"); j_theta)
        theta = j_theta->get<double>()/180*std::numbers::pi;

      solver.load(*j("solver"));
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
