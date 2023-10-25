#pragma once

#include "pore_network_image.hpp"

namespace xpm {
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

  class phase_state
  {
    dpl::strong_vector<net_t, phase_config> config_;
    std::vector<phase_config> config_t_;

    dpl::strong_vector<net_t, double> local_;
    std::vector<double> local_t_;

    static inline constexpr auto r_cap_mobile_ = std::numeric_limits<double>::min();

  public:
    double r_cap_global;

    void resize(net_t mv, std::size_t t) {
      config_.resize(mv);
      config_t_.resize(t);

      local_.resize(mv, r_cap_mobile_);
      local_t_.resize(t, r_cap_mobile_);
    }

    auto& config(net_t i) const { return config_[i]; }
    auto& config(net_t i) { return config_[i]; }

    auto& config(std::size_t i) const { return config_t_[i]; }
    auto& config(std::size_t i) { return config_t_[i]; }

    auto& local(net_t i) const { return local_[i]; }
    auto& local(net_t i) { return local_[i]; }

    auto& local(std::size_t i) const { return local_t_[i]; }
    auto& local(std::size_t i) { return local_t_[i]; }

    bool mobile(const auto i) const {
      return local(i) == r_cap_mobile_;  // NOLINT(clang-diagnostic-float-equal)
    }

    auto r_cap(const auto i) const {
      auto l = local(i);
      return l == r_cap_mobile_ ? r_cap_global : l;  // NOLINT(clang-diagnostic-float-equal)
    }
  };

  class invasion_task {
    template <bool phase1>
    class filter_functor
    {
    public:
      const pore_network_image* pni;
      const phase_state* state;
      bool darcy_filter;

      bool operator()(macro_t i) const {
        if constexpr (phase1)
          return state->config(pni->net(i)).phase() == phase_config::phase1();
        else
          return true /*state->macro_voxel[pni->net(i)].phase() == phase_config::phase0()*/; // TODO
      }

      bool operator()(std::size_t i) const {
        if constexpr (phase1)
          return state->config(i).phase() == phase_config::phase1();
        else
          return true /*state->throat[i].phase() == phase_config::phase0()*/; // TODO
      }

      bool operator()(voxel_t i) const {
        return darcy_filter;

        // if constexpr (phase1)
        //   return state->macro_voxel[pni->net(i)].phase() == phase_config::phase1();
        // else
        //   return true /*state->macro_voxel[pni->net(i)].phase() == phase_config::phase0()*/;
      }
    };


    class phase1_coef_functor
    {
    public:
      const pore_network_image* pni;
      const phase_state* state;
      double theta;
      double darcy_perm;

      auto operator()(macro_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;
        auto area = eq_tr::area(attrib::r_ins(pni->pn(), i));
        return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(pni->net(i)))/area);
      }

      auto operator()(size_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;
        auto area = eq_tr::area(attrib::r_ins(pni->pn(), i));
        return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(i))/area);
      }

      auto operator()(voxel_t) const {
        return darcy_perm; // TODO
      }
    };


    class phase0_coef_functor
    {
    public:
      const pore_network_image* pni;
      const phase_state* state;
      double theta;
      double darcy_perm;

      double operator()(std::size_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;

        return
          state->config(i).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(i)))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni->pn(), i)));
      }

      double operator()(macro_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;

        return
          state->config(pni->net(i)).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(pni->net(i))))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni->pn(), i)));
      }

      double operator()(voxel_t) const {
        return darcy_perm; // TODO
      }
    };

    using vertex_t = dpl::graph::dc_graph::vertex_t;
    using edge_t = dpl::graph::dc_graph::edge_t;
    using et_ptr = dpl::graph::et_traits::node_ptr;

    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;
    const startup_settings* settings_;

    double total_pore_volume_;

    dpl::graph::dc_graph g_;
    std::unordered_map<edge_t, std::size_t> de_to_throat_idx_;
    std::unique_ptr<edge_t[]> throat_idx_to_de_;

    dpl::graph::dc_context<dpl::graph::dc_traits> context_;

    phase_state state_;
    std::vector<double> kr_recorded_;


    dpl::vector2d last_pc_point_;
    double pc_max_step_;

    struct curves
    {
      std::vector<dpl::vector2d> pc;
      std::vector<dpl::vector3d> kr;

      bool add_pc_point(const dpl::vector2d& p) {
        if (pc.size() > 1
          && std::abs(pc[pc.size() - 1].y() - p.y()) < 1e-6
          && std::abs(pc[pc.size() - 2].y() - p.y()) < 1e-6) {
          pc[pc.size() - 1].x() = p.x();
          return false;
        }

        pc.push_back(p);
        return true;
      }
    } primary_,
      secondary_;

    idx1d_t progress_idx_ = 0;

    bool darcy_invaded_ = false;
    bool finished_ = false;

    

    

  public:
    explicit invasion_task(pore_network_image& pni, const startup_settings& settings) :
      pni_{&pni},
      pn_{&pni.pn()},
      img_{&pni.img()},
      settings_{&settings} {}

    auto& graph() {
      return g_;
    }

    auto& state() {
      return state_;
    }

    auto& primary() {
      return primary_;
    }

    auto& secondary() {
      return secondary_;
    }

    auto progress_idx() const {
      return progress_idx_;
    }

    bool finished() const {
      return finished_;
    }

    void init() {
      state_.resize(pni_->connected_count(), pn_->throat_count());

      primary_.pc.reserve(500);
      primary_.kr.reserve(500);
      secondary_.pc.reserve(500);
      secondary_.kr.reserve(500);

      total_pore_volume_ = pni_->total_pore_volume(settings_->darcy_poro);


      {
        using namespace std::chrono;
        std::cout << "graph...";
        auto t0 = high_resolution_clock::now();
        std::tie(g_, de_to_throat_idx_, throat_idx_to_de_) = pni_->generate_dc_graph<true>();
        std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
          duration_cast<seconds>(high_resolution_clock::now() - t0).count(), g_.vertex_count(), g_.edge_count());
      }

      std::cout << "full decremental connectivity... (async)\n\n";
    }

    void generate_euler_tour() {
      using namespace std::chrono;
      using clock = high_resolution_clock;

      {
        std::cout << "euler tour...";
        auto t0 = clock::now();
        context_.init_with_dfs(g_, dpl::graph::dc_traits{g_});
        std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
      }
    }

    double calc_relative(const startup_settings::input_curves& curve, std::span<const dpl::vector2d> pc_inv, double theta, auto phase1) {
      auto darcy_sw = curve.kr0.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : solve(pc_inv, 1/state_.r_cap_global, dpl::extrapolant::flat);

      double kr = 0.0;

      if constexpr (phase1) {
        if (!curve.kr1.empty())
          kr = std::max(0.0, solve(std::span{curve.kr1}, darcy_sw, dpl::extrapolant::flat));
      }
      else {
        if (!curve.kr0.empty())
          kr = std::max(0.0, solve(std::span{curve.kr0}, darcy_sw, dpl::extrapolant::flat));
      }

      static constexpr auto kr_threshold = 1e-3;

      bool darcy_filter = kr > kr_threshold;

      filter_functor<phase1> filter{pni_, &state_, darcy_filter};
      auto [nrows, mapping] = pni_->generate_mapping(*settings_->solver.decomposition, filter);

      std::conditional_t<phase1, phase1_coef_functor, phase0_coef_functor> term{pni_, &state_, theta, kr*settings_->darcy_perm};

      auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);

      using namespace std::chrono;

      auto hash = std::hash<std::string_view>{}(
        std::string_view{reinterpret_cast<char*>(input.values.get()), nvalues*sizeof(HYPRE_Complex)});

      std::cout << fmt::format("ph{} | K: {:.2e} mD | hash: {:x}", // kr: {:.2e}, dar_sw: {:.4f}, 
        int{phase1}, kr*settings_->darcy_perm/presets::darcy_to_m2, hash); /*kr_val, darcy_sw, */

      auto found = pressure_cache::cache().find(hash);

      if (found == pressure_cache::cache().end()) {
        dpl::strong_vector<net_t, double> pressure(pni_->connected_count());
        auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);

        {
          // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);

          auto t0 = high_resolution_clock::now();
          dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_->solver.tolerance, settings_->solver.max_iterations);
          auto t1 = high_resolution_clock::now();
          std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings_->solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
          auto t2 = high_resolution_clock::now();
          std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);

          std::cout << fmt::format(" | save: {} s, solve: {} s\n",
            duration_cast<seconds>(t1 - t0).count(),
            duration_cast<seconds>(t2 - t1).count());
        }
        
        for (HYPRE_BigInt i = 0; i < nrows; ++i)
          pressure[mapping.backward[i]] = decomposed_pressure[i];

        auto [inlet, outlet] = pni_->flow_rates(pressure, term, filter);

        pressure_cache::cache()[hash] = inlet;
        pressure_cache::save();

        return inlet;
      }

      std::cout << '\n';
      return found->second;
    }

    void launch_secondary(double absolute_rate, double theta) {
      using eq_tr = hydraulic_properties::equilateral_triangle;
      using namespace attrib;


      auto pc_inv = settings_->secondary.pc;
      std::ranges::reverse(pc_inv);
      for (auto& p : pc_inv)
        std::swap(p.x(), p.y());
      pc_inv.resize(std::ranges::unique(pc_inv, {}, [](const dpl::vector2d& p) { return p.x(); }).begin() - pc_inv.begin());

      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim()).prod()*settings_->darcy_poro;


      double darcy_r_cap = std::numeric_limits<double>::quiet_NaN();
      if (!pc_inv.empty())
        darcy_r_cap = 1/settings_->secondary.pc.back().y();


      dpl::strong_vector<net_t, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue<false> queue;

      auto area_corner_mult = 0.0;

      for (macro_t i{0}; i < pn_->node_count(); ++i)
        if (pni_->connected(i) && state_.config(pni_->net(i)).layout() == phase_config::bulk_films()) {
          area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));
          queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
        }
      
      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = adj(pn_, i); pni_->connected(l) && state_.config(i).layout() == phase_config::bulk_films()) {
          area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));

          if (pn_->inlet() == r) {
            queue.insert(i, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, i)));
            explored_throat[i] = true;
          }
          else
            queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
        }

      idx1d_t darcy_count = 0;


      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_t idx1d{0};

        for (k = 0; k < img_->dim().z(); ++k)
          for (j = 0; j < img_->dim().y(); ++j)
            for (i = 0; i < img_->dim().x(); ++i, ++idx1d)
              if (pni_->connected(idx1d)/* && filter(idx1d)*/) { // TODO
                ++darcy_count;
                queue.insert(idx1d, darcy_r_cap);

                // if (auto velem = img_->velem[idx1d]; velem/* && filter(velem)*/) // macro-darcy
                //   builder.reserve(velem, idx1d);
                //
                // dpl::sfor<3>([&](auto d) {
                //   if (ijk[d] < img_->dim()[d] - 1)
                //     if (auto adj_idx1d = idx1d + img_->idx_map(d); img_->phase[adj_idx1d] == presets::microporous && filter(adj_idx1d)) // darcy-darcy
                //       builder.reserve(idx1d, adj_idx1d);
                // });
              }
      }

      

      auto inv_volume_coef0 = 0.0;

      auto eval_inv_volume = [&] {
        auto macro = area_corner_mult*eq_tr::area_corners(theta, state_.r_cap_global) + inv_volume_coef0;
        if (!pc_inv.empty())
          macro += solve(std::span<const dpl::vector2d>{pc_inv}, 1/state_.r_cap_global, dpl::extrapolant::flat)*unit_darcy_pore_volume*darcy_count;
        return macro;

        // return area_corner_mult*eq_tr::area_corners(theta, state_.r_cap_global) + inv_volume_coef0;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };

      using namespace dpl::graph;

      net_t outlet_idx = pni_->connected_count();
      dc_traits traits{g_};
      et_ptr outlet_entry = traits.get_entry(vertex_t(*outlet_idx)); 

      auto freeze_cluster = [&](et_ptr hdr) {
        auto area_corners = eq_tr::area_corners(theta, state_.r_cap_global);

        for (et_ptr et : range<et_traits>(hdr))
          if (dc_traits::is_loop_edge(et)) { // TODO: CHECK IF MACRO/DARCY      if (pni_->is_macro(b_net_idx)) 
            vertex_t v = dc_traits::get_vertex(et);

            traits.set_entry(v, nullptr);

            if (net_t net_idx{*v}; state_.config(net_idx).layout() == phase_config::bulk_films()) {  // NOLINT(clang-diagnostic-float-equal)
              state_.local(net_idx) = state_.r_cap_global;

              if (pni_->is_macro(net_idx)) {
                auto macro_idx = pni_->macro(net_idx);

                auto mult = volume(pn_, macro_idx)/eq_tr::area(r_ins(pn_, macro_idx));
                area_corner_mult -= mult;
                inv_volume_coef0 += area_corners*mult;  
              }
              else {
                --darcy_count;
                inv_volume_coef0 += solve(std::span<const dpl::vector2d>{pc_inv}, 1/state_.r_cap_global, dpl::extrapolant::flat)*unit_darcy_pore_volume;  
              }
            }

            for (edge_t vu : g_.edges(v)) {
              auto net_idx = net_t{*target(vu, g_)};

              if (pni_->is_macro(net_idx)) {
                auto t_idx = de_to_throat_idx_[vu]; // TODO: CHECK edge is MACRO-MACRO throat

                if (state_.config(t_idx).layout() == phase_config::bulk_films() && state_.mobile(t_idx)) {

                  traits.set_null_entry(vu);
                  traits.set_null_entry(opposite(vu, g_));  

                  state_.local(t_idx) = state_.r_cap_global;

                  auto mult = volume(pn_, t_idx)/eq_tr::area(r_ins(pn_, t_idx));
                  area_corner_mult -= mult;
                  inv_volume_coef0 += area_corners*mult;
                }
              } 
              else {
                traits.set_null_entry(vu);
                traits.set_null_entry(opposite(vu, g_));  // TODO: REMOVE - REDUNDANT
              }
            }
          }
      };



      double last_kr_sw = 0.0;
      secondary_.kr.emplace_back(0, 0, 1);
      // secondary_.add_pc_point(last_pc_point_); // TODO





      // auto rel_calc = [=, this, &pc_inv](auto phase1) {
      //   auto darcy_sw = settings_->secondary.kr0.empty()
      //     ? std::numeric_limits<double>::quiet_NaN()
      //     : solve(std::span<const dpl::vector2d>{pc_inv}, 1/state_.r_cap_global, dpl::extrapolant::flat);
      //
      //   double kr = 0.0;
      //
      //   if constexpr (phase1) {
      //     if (!settings_->secondary.kr1.empty())
      //       kr = std::max(0.0, solve(std::span{settings_->secondary.kr1}, darcy_sw, dpl::extrapolant::flat));
      //   }
      //   else {
      //     if (!settings_->secondary.kr0.empty())
      //       kr = std::max(0.0, solve(std::span{settings_->secondary.kr0}, darcy_sw, dpl::extrapolant::flat));
      //   }
      //
      //   static constexpr auto kr_threshold = 1e-3;
      //
      //   bool darcy_filter = kr > kr_threshold;
      //
      //   filter_functor<phase1> filter{pni_, &state_, darcy_filter};
      //   auto [nrows, mapping] = pni_->generate_mapping(*settings_->solver.decomposition, filter);
      //
      //   std::conditional_t<phase1, phase1_coef_functor, phase0_coef_functor> term{pni_, &state_, theta, kr*settings_->darcy_perm};
      //
      //   auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);
      //
      //   using namespace std::chrono;
      //
      //   auto hash = std::hash<std::string_view>{}(
      //     std::string_view{reinterpret_cast<char*>(input.values.get()), nvalues*sizeof(HYPRE_Complex)});
      //
      //   std::cout << fmt::format("ph{} | K: {:.2e} mD | hash: {:x}", // kr: {:.2e}, dar_sw: {:.4f}, 
      //     int{phase1}, kr*settings_->darcy_perm/presets::darcy_to_m2, hash); /*kr_val, darcy_sw, */
      //
      //   auto found = pressure_cache::cache().find(hash);
      //
      //   if (found == pressure_cache::cache().end()) {
      //     dpl::strong_vector<net_t, double> pressure(pni_->connected_count());
      //     auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);
      //
      //     {
      //       // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);
      //
      //       auto t0 = high_resolution_clock::now();
      //       dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_->solver.tolerance, settings_->solver.max_iterations);
      //       auto t1 = high_resolution_clock::now();
      //       std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings_->solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
      //       auto t2 = high_resolution_clock::now();
      //       std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);
      //
      //       std::cout << fmt::format(" | save: {} s, solve: {} s\n",
      //         duration_cast<seconds>(t1 - t0).count(),
      //         duration_cast<seconds>(t2 - t1).count());
      //     }
      //     
      //     for (HYPRE_BigInt i = 0; i < nrows; ++i)
      //       pressure[mapping.backward[i]] = decomposed_pressure[i];
      //
      //     auto [inlet, outlet] = pni_->flow_rates(pressure, term, filter);
      //
      //     pressure_cache::cache()[hash] = inlet;
      //     pressure_cache::save();
      //
      //     return inlet;
      //   }
      //
      //   std::cout << '\n';
      //   return found->second;
      // };


      auto rel_calc_report = [&](double sw) {
        auto phase1_connected = [&] {
          for (auto [l, r] : pn_->throat_.span(adj))
            if (r == pn_->inlet() && pni_->connected(l)) // macro-inlet
              if (traits.get_entry(vertex_t(*pni_->net(l))))
                return true;

          {
            idx3d_t ijk;
            auto& [i, j, k] = ijk;
          
            for (k = 0; k < img_->dim().z(); ++k)
              for (j = 0; j < img_->dim().y(); ++j)
                if (voxel_t idx1d{img_->idx_map(0, j, k)}; pni_->connected(idx1d)/*img_->phase[idx1d] == presets::microporous*/) // darcy-inlet
                  if (traits.get_entry(vertex_t(*pni_->net(idx1d))))
                    return true;
          }

          return false;
        };


        
        auto def = calc_relative(settings_->secondary, pc_inv, theta, std::false_type{});//rel_calc(std::false_type{});
        auto inv = phase1_connected() ? calc_relative(settings_->secondary, pc_inv, theta, std::true_type{})/*rel_calc(std::true_type{})*/ : 0;
      
        secondary_.kr.emplace_back(sw, def/absolute_rate, inv/absolute_rate);
      };





      if (!pc_inv.empty()) {
        auto max_pc = primary_.pc.back().y();

        state_.r_cap_global = 1/max_pc;
        secondary_.add_pc_point(eval_pc_point());
        
        {
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_pc) - std::log10(1/darcy_r_cap))/steps);

          for (auto i = 0; i < steps; ++i) {
            state_.r_cap_global *= step;
            last_pc_point_ = eval_pc_point();
            secondary_.add_pc_point(last_pc_point_);
          }
        }

        ++progress_idx_;

        if (!settings_->secondary.kr0.empty()) {
          auto steps = 6;

          for (auto i = 1; i < steps; ++i) {
            state_.r_cap_global = 1/solve(std::span{settings_->secondary.pc}, 1.*i/steps, dpl::extrapolant::flat);
            auto sw = eval_inv_volume()/total_pore_volume_;
            last_kr_sw = sw;
            rel_calc_report(sw);
            // rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
            ++progress_idx_;
          }
        }
      }


      auto inv_idx_ = 0;

      for (; !queue.empty(); ++progress_idx_) {
        // if (++inv_idx_ < 5000)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{1});



        if (auto sw = eval_inv_volume()/total_pore_volume_; (sw - last_kr_sw) > 0.075) {
          last_kr_sw = sw;
          rel_calc_report(sw);
        }

        if (auto pc_point = eval_pc_point();
          std::abs(last_pc_point_.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point_.y() - pc_point.y()) > pc_max_step_) {
          last_pc_point_ = pc_point;
          secondary_.add_pc_point(pc_point);
        }

        auto [elem, local_idx, r_cap] = queue.front();

        

        queue.pop();

        if (elem == displ_elem::macro) {
          macro_t macro_idx(local_idx);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net_idx = pni_->net(macro_idx);

          if (state_.config(net_idx).phase() == phase_config::phase0())
            continue;

          if (traits.get_entry(vertex_t(*net_idx))) {
            state_.r_cap_global = std::max(r_cap, state_.r_cap_global);
            state_.config(net_idx) = phase_config{phase_config::phase0()};
            
            area_corner_mult -= volume(pn_, macro_idx)/eq_tr::area(r_ins(pn_, macro_idx));
            inv_volume_coef0 += volume(pn_, macro_idx);

            {
              vertex_t v(*net_idx);

              for (edge_t vu : g_.edges(v))
                if (!traits.is_null_entry(vu) && !traits.is_tree_edge(vu))
                  context_.non_tree_edge_remove(vu);

              for (edge_t vu : g_.edges(v)) {
                if (net_t u_net_idx{*target(vu, g_)}; u_net_idx == outlet_idx || pni_->is_macro(u_net_idx)) { // macro-macro
                  if (auto t_idx = de_to_throat_idx_[vu]; !explored_throat[t_idx]) {
                    queue.insert(t_idx, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, t_idx)));
                    explored_throat[t_idx] = true;
                  }
                }
                else { // macro-darcy
                  // TODO Darcy 
                }

                if (!traits.is_null_entry(vu))
                  if (!context_.tree_edge_split_and_reconnect(vu))
                    if (et_ptr hdr = et_algo::get_header(traits.get_entry(target(vu, g_))); hdr != et_algo::get_header(outlet_entry))
                      freeze_cluster(hdr);
              }

              traits.set_entry(v, nullptr);
            }
          }
        }
        else if (elem == displ_elem::throat) {
          if (state_.config(local_idx).phase() == phase_config::phase0())
            continue;

          auto [l, r] = adj(pn_, local_idx);

          auto l_net = pni_->net(l);

          et_ptr l_entry = traits.get_entry(vertex_t(*l_net));
          et_ptr r_entry =
            r == pn_->inlet() ? nullptr :
            r == pn_->outlet() ? outlet_entry :
            traits.get_entry(vertex_t(*pni_->net(r)));

          if (l_entry || r_entry) {
            if (l_entry && !explored[l_net]) {
              queue.insert(l, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, l)));
              explored[l_net] = true;
            }

            if (pn_->inner_node(r) && r_entry)
              if (auto r_net = pni_->net(r); !explored[r_net]) {
                queue.insert(r, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, r)));
                explored[r_net] = true;
              }

            state_.config(local_idx) = phase_config{phase_config::phase0()};
            state_.r_cap_global = std::max(r_cap, state_.r_cap_global);
            
            area_corner_mult -= volume(pn_, local_idx)/eq_tr::area(r_ins(pn_, local_idx));
            inv_volume_coef0 += volume(pn_, local_idx);

            if (edge_t de = throat_idx_to_de_[local_idx]; !traits.is_null_entry(de))
              if (traits.is_tree_edge(de)) {
                if (!context_.tree_edge_split_and_reconnect(de)) {
                  et_ptr hdr = et_algo::get_header(traits.get_entry(target(de, g_)));
                  if (hdr == et_algo::get_header(outlet_entry))
                    hdr = et_algo::get_header(traits.get_entry(target(opposite(de, g_), g_)));

                  freeze_cluster(hdr);
                }
              }
              else
                context_.non_tree_edge_remove(de);
          }
        }
        else { // voxel
          voxel_t voxel_idx(local_idx);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net_idx = pni_->net(voxel_idx);
          
          // if (traits.get_entry(vertex_t(*net_idx))) 
          {
            state_.r_cap_global = std::max(r_cap, state_.r_cap_global);
            state_.config(net_idx) = phase_config{phase_config::phase0()};

            vertex_t v(*net_idx);
              
            for (edge_t vu : g_.edges(v))
              if (!traits.is_null_entry(vu) && !traits.is_tree_edge(vu))
                context_.non_tree_edge_remove(vu);
              
            for (edge_t vu : g_.edges(v))
              if (!traits.is_null_entry(vu) && !context_.tree_edge_split_and_reconnect(vu))
                if (et_ptr hdr = et_algo::get_header(traits.get_entry(target(vu, g_))); hdr != et_algo::get_header(outlet_entry))
                  freeze_cluster(hdr);

            traits.set_entry(v, nullptr);
          }
        }
      }


      secondary_.add_pc_point(eval_pc_point());
      rel_calc_report(eval_inv_volume()/total_pore_volume_);
      ++progress_idx_;
    }

    void launch_primary(
      double absolute_rate,
      double theta,
      std::span<const dpl::vector2d> pc_inv) {

      pressure_cache::load();

      using eq_tr = hydraulic_properties::equilateral_triangle;
      using namespace attrib;

      
      double darcy_r_cap;

      if (pc_inv.empty()) {
        auto min_r_cap_throat = std::numeric_limits<double>::max();

        for (std::size_t i{0}; i < pn_->throat_count(); ++i)
          if (auto [l, r] = adj(pn_, i); pn_->inner_node(r) && pni_->connected(l))
            min_r_cap_throat = std::min(min_r_cap_throat, r_ins(pn_, i));
        min_r_cap_throat = 0.95*eq_tr::r_cap_piston_with_films(theta, min_r_cap_throat);
      
        darcy_r_cap = min_r_cap_throat;
      }
      else
        darcy_r_cap = 1/pc_inv.front().x();

      pc_max_step_ = 0.075/darcy_r_cap;

      dpl::strong_vector<net_t, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue<true> queue;

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = adj(pn_, i); pn_->inlet() == r) // inlet macro-macro TODO: voxels inlet needed
          if (pni_->connected(l)) {
            queue.insert(i, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, i)));
            explored_throat[i] = true;
          }

      using namespace dpl::graph;

      net_t outlet_idx = pni_->connected_count();

      dpl::vector3d inv_volume_coefs{0};
      idx1d_t inv_darcy_count{0};

      state_.r_cap_global = queue.front().radius_cap;

      last_pc_point_ = {2, -1};
      double last_kr_sw = 1.0;
      primary_.kr.emplace_back(1, 1, 0);

      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim()).prod()*settings_->darcy_poro;

      auto eval_inv_volume = [&] {
        auto macro = inv_volume_coefs[0] + inv_volume_coefs[2]*state_.r_cap_global*state_.r_cap_global;
        if (!pc_inv.empty())
          macro += (1 - solve(pc_inv, 1/state_.r_cap_global, dpl::extrapolant::flat))*unit_darcy_pore_volume*inv_darcy_count;
        return macro;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{1 - eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };
      
      auto rel_calc_report = [=, this](double sw) {
        auto ph0 = calc_relative(settings_->primary, pc_inv, theta, std::false_type{});
        auto ph1 = calc_relative(settings_->primary, pc_inv, theta, std::true_type{}); // TODO: calculate when breakthrough
        
        primary_.kr.emplace_back(sw, ph0/absolute_rate, ph1/absolute_rate);
      };

      for (progress_idx_ = 0; !queue.empty(); ++progress_idx_) {
        // if (inv_idx_ < 35)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{150});
        // else if (inv_idx_ < 100)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{50});
        // else if (inv_idx_ < 1000)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{25});

        if (auto sw = (1 - eval_inv_volume()/total_pore_volume_); (last_kr_sw - sw) > 0.075) {
          last_kr_sw = sw;
          rel_calc_report(sw);
        }

        if (auto pc_point = eval_pc_point();
          std::abs(last_pc_point_.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point_.y() - pc_point.y()) > pc_max_step_) {
          last_pc_point_ = pc_point;
          primary_.add_pc_point(pc_point);
        }

        auto [elem, local_idx, r_cap] = queue.front();

        queue.pop();

        if (elem == displ_elem::macro) {
          auto macro_idx = macro_t(local_idx);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net_idx = pni_->net(macro_idx);

          state_.config(net_idx) = phase_config::phase1_bulk_films();
          state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

          inv_volume_coefs[0] += volume(pn_, macro_idx);
          inv_volume_coefs[2] -= volume(pn_, macro_idx)*eq_tr::area_films(theta)/eq_tr::area(r_ins(pn_, macro_idx));
          
          for (edge_t vu : g_.edges(vertex_t(*net_idx)))
            if (net_t u_net_idx{*target(vu, g_)}; u_net_idx == outlet_idx || pni_->is_macro(u_net_idx)) { // macro-macro
              if (auto t_idx = de_to_throat_idx_[vu]; !explored_throat[t_idx]) {
                queue.insert(t_idx, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, t_idx)));
                explored_throat[t_idx] = true;
              }
            }
            else { // macro-darcy
              if (!explored[u_net_idx]) {
                auto b_voxel_idx = pni_->voxel(u_net_idx);
                queue.insert(b_voxel_idx, darcy_r_cap);
                explored[u_net_idx] = true;
              }
            }
        }
        else if (elem == displ_elem::voxel) {
          auto net_idx = pni_->net(voxel_t(local_idx));

          darcy_invaded_ = true;
          ++inv_darcy_count;

          state_.config(net_idx) = phase_config::phase1_bulk_films(); // TODO
          state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

          for (auto vu : g_.edges(vertex_t(*net_idx))) {
            if (net_t u_net_idx{*target(vu, g_)}; u_net_idx != outlet_idx && !explored[u_net_idx]) {
              if (pni_->is_macro(u_net_idx)) { // macro
                auto v_macro_idx = pni_->macro(u_net_idx);
                queue.insert(v_macro_idx, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, v_macro_idx)));
                explored[u_net_idx] = true;
              }
              else { // darcy
                queue.insert(pni_->voxel(u_net_idx), darcy_r_cap);
                explored[u_net_idx] = true;
              }
            }
          }
        }
        else { // throat
          auto [l, r] = adj(pn_, local_idx);

          state_.config(local_idx) = phase_config::phase1_bulk_films();
          state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

          inv_volume_coefs[0] += volume(pn_, local_idx);
          inv_volume_coefs[2] -= volume(pn_, local_idx)*eq_tr::area_films(theta)/eq_tr::area(r_ins(pn_, local_idx));

          if (pn_->inner_node(r))
            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(r, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, r)));
              explored[r_net] = true;
            }

          if (auto l_net = pni_->net(l); !explored[l_net]) {
            queue.insert(l, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, l)));
            explored[l_net] = true;
          }
        }
      }


      if (!pc_inv.empty()) {
        auto end_pc = 1/state_.r_cap_global;
        auto max_darcy_pc = pc_inv.back().x();

        primary_.add_pc_point(eval_pc_point());

        { /* Capillary pressure */
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_darcy_pc) - std::log10(end_pc))/steps);

          for (auto i = 0; i < steps; ++i) {
            state_.r_cap_global /= step;
            primary_.add_pc_point(eval_pc_point());
          }
        }

        ++progress_idx_;

        state_.r_cap_global = 1/end_pc;

        if (!settings_->primary.kr0.empty()) {
          auto steps = 6;

          for (auto i = 1; i < steps; ++i) {
            state_.r_cap_global = 1/solve(std::span{settings_->primary.pc}, 1 - 1.*i/steps, dpl::extrapolant::flat);
            rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
            ++progress_idx_;
          }
        }
      }

      primary_.kr.emplace_back(0, 0, 1);
      ++progress_idx_;

      finished_ = true;

      {
        std::cout << '\n';
        generate_euler_tour();
        launch_secondary(absolute_rate, theta);
      }
    }
  };
}