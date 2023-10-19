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
    dpl::strong_vector<net_tag, phase_config> config_;
    std::vector<phase_config> config_t_;

    dpl::strong_vector<net_tag, double> local_;
    std::vector<double> local_t_;

    static inline constexpr auto r_cap_mobile_ = std::numeric_limits<double>::min();

  public:
    double r_cap_global;

    void resize(net_idx_t mv, std::size_t t) {
      config_.resize(mv);
      config_t_.resize(t);

      local_.resize(mv, r_cap_mobile_);
      local_t_.resize(t, r_cap_mobile_);
    }

    auto& config(net_idx_t i) const { return config_[i]; }
    auto& config(net_idx_t i) { return config_[i]; }

    auto& config(std::size_t i) const { return config_t_[i]; }
    auto& config(std::size_t i) { return config_t_[i]; }

    auto& local(net_idx_t i) const { return local_[i]; }
    auto& local(net_idx_t i) { return local_[i]; }

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

      bool operator()(macro_idx_t i) const {
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

      bool operator()(voxel_idx_t i) const {
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

      auto operator()(macro_idx_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle_properties;
        auto area = eq_tr::area(attribs::r_ins(pni->pn(), i));
        return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(pni->net(i)))/area);
      }

      auto operator()(size_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle_properties;
        auto area = eq_tr::area(attribs::r_ins(pni->pn(), i));
        return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(i))/area);
      }

      auto operator()(voxel_idx_t) const {
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
        using eq_tr = hydraulic_properties::equilateral_triangle_properties;

        return
          state->config(i).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(i)))
            : eq_tr::conductance_single(eq_tr::area(attribs::r_ins(pni->pn(), i)));
      }

      double operator()(macro_idx_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle_properties;

        return
          state->config(pni->net(i)).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(pni->net(i))))
            : eq_tr::conductance_single(eq_tr::area(attribs::r_ins(pni->pn(), i)));
      }

      double operator()(voxel_idx_t) const {
        return darcy_perm; // TODO
      }
    };






    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;
    const startup_settings* settings_;

    double total_pore_volume_;

    dpl::graph::dc_graph graph_;
    std::unordered_map<dpl::graph::dc_graph::edge_descriptor, std::size_t> de_to_throat_;
    std::unique_ptr<dpl::graph::dc_graph::edge_descriptor[]> throat_to_de_;

    dpl::graph::dc_context<dpl::graph::dc_properties> dc_context_;

    phase_state state_;

    dpl::vector2d last_pc_point_;
    double pc_max_step_;
    std::vector<dpl::vector2d> pc_curve_;
    std::vector<dpl::vector3d> kr_curves_;

    idx1d_t progress_idx_ = 0;

    bool darcy_invaded_ = false;
    bool finished_ = false;

    void add_to_pc_curve(const dpl::vector2d& p) {
      if (pc_curve_.size() > 1
        && std::abs(pc_curve_[pc_curve_.size() - 1].y() - p.y()) < 1e-6
        && std::abs(pc_curve_[pc_curve_.size() - 2].y() - p.y()) < 1e-6)
        pc_curve_[pc_curve_.size() - 1].x() = p.x();
      else
        pc_curve_.push_back(p);
    }

  public:
    explicit invasion_task(pore_network_image& pni, const startup_settings& settings) :
      pni_{&pni},
      pn_{&pni.pn()},
      img_{&pni.img()},
      settings_{&settings} {}

    auto& graph() {
      return graph_;
    }

    auto& state() {
      return state_;
    }

    auto& pc_curve() const {
      return pc_curve_;
    }

    auto& kr_curves() const {
      return kr_curves_;
    }

    auto progress_idx() const {
      return progress_idx_;
    }

    bool finished() const {
      return finished_;
    }

    void init() {
      state_.resize(pni_->connected_count(), pn_->throat_count());

      pc_curve_.reserve(1000);
      kr_curves_.reserve(1000);

      total_pore_volume_ = pni_->total_pore_volume(settings_->darcy_poro);

      std::cout << "full decremental connectivity... (async)\n";
    }

    void generate_graph() {
      using namespace std::chrono;
      using clock = high_resolution_clock;

      {
        std::cout << "graph...";
        auto t0 = clock::now();
        std::tie(graph_, de_to_throat_, throat_to_de_) = pni_->generate_dc_graph<true>();
        std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
          duration_cast<seconds>(clock::now() - t0).count(), graph_.vertex_count(), graph_.edge_count());
      }

      {
        std::cout << "euler tour...";
        auto t0 = clock::now();
        dc_context_.init_with_dfs(graph_, dpl::graph::dc_properties{graph_});
        std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
      }
    }

    void launch_secondary(
      double absolute_rate,
      double theta_rec,
      double theta_adv) {

      displ_queue<false> queue;

      using eq_tr = hydraulic_properties::equilateral_triangle_properties;
      using namespace attribs;

      auto area_corner_mult = 0.0;

      for (macro_idx_t i{0}; i < pn_->node_count(); ++i) {
        if (pni_->connected(i) && state_.config(pni_->net(i)).layout() == phase_config::bulk_films()) {
          queue.insert(displ_elem::macro, *i, eq_tr::r_cap_snap_off(theta_adv, r_ins(pn_, i)));
          area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));
        }
      }

      for (std::size_t i{0}; i < pn_->throat_count(); ++i) {
        auto [l, r] = pn_->throat_[adj][i];
        if (pni_->connected(l) && pn_->inner_node(r) && state_.config(i).layout() == phase_config::bulk_films()) {
          queue.insert(displ_elem::throat, i, eq_tr::r_cap_snap_off(theta_adv, r_ins(pn_, i)));
          area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));
        }
      }

      auto inv_volume_coef0 = 0.0;

      auto eval_area_corners = [&] {
        return eq_tr::area_corners(theta_adv, state_.r_cap_global);
      };

      auto eval_inv_volume = [&] {
        return area_corner_mult*eval_area_corners() + inv_volume_coef0;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };

      using namespace dpl::graph;

      auto freeze_cluster = [&](et_traits::node_ptr hdr) {
        auto area_corners = eval_area_corners();

        for (auto et : range<et_traits>(hdr))
          if (dc_properties::is_loop_edge(et)) { // TODO: CHECK IF MACRO/DARCY      if (pni_->is_macro(b_net_idx)) 
            auto a = dc_properties::get_vertex(et);

            net_idx_t net_idx{a};
            
            if (state_.mobile(net_idx)) {  // NOLINT(clang-diagnostic-float-equal)
              state_.local(net_idx) = state_.r_cap_global;

              auto inner_macro_idx = pni_->macro(net_idx);

              auto mult = volume(pn_, inner_macro_idx)/eq_tr::area(r_ins(pn_, inner_macro_idx));
              area_corner_mult -= mult;
              inv_volume_coef0 += area_corners*mult;
            }

            for (auto ab : graph_.edges(a)) {
              auto b = target(ab, graph_);
              auto t_idx = de_to_throat_[ab]; // TODO: CHECK edge is MACRO-MACRO throat
              
              if (pni_->is_macro(net_idx_t{b}) && state_.config(t_idx).layout() == phase_config::bulk_films()) {
                if (state_.mobile(t_idx)) { // NOLINT(clang-diagnostic-float-equal)
                  state_.local(t_idx) = state_.r_cap_global;

                  auto mult = volume(pn_, t_idx)/eq_tr::area(r_ins(pn_, t_idx));
                  area_corner_mult -= mult;
                  inv_volume_coef0 += area_corners*mult;
                }
              }
            }
          }
      };

      // pn_->node_[volume][local_idx]*
      //       -eq_tr::area_of_films(theta)/eq_tr::area(pn_->node_[r_ins][local_idx])


      net_idx_t outlet_idx{*pni_->connected_count()};
      dc_properties dc_props{graph_};
      et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);

      double last_kr_sw = 0.0;






      auto rel_calc = [=, this](auto phase1) {
        // auto darcy_sw = settings_->darcy_kr0.empty()
        //   ? std::numeric_limits<double>::quiet_NaN()
        //   : solve(darcy_pc_inv, 1/state_.r_cap, dpl::extrapolant::flat);
      
        // double kr_val = 0.0;
        //
        // if constexpr (phase1) {
        //   if (!settings_->darcy_kr1.empty())
        //     kr_val = std::max(0.0, solve(std::span{settings_->darcy_kr1}, darcy_sw, dpl::extrapolant::flat));
        // }
        // else {
        //   if (!settings_->darcy_kr0.empty())
        //     kr_val = std::max(0.0, solve(std::span{settings_->darcy_kr0}, darcy_sw, dpl::extrapolant::flat));
        // }
      
        // static constexpr auto perm_threshold = 1e-6*presets::darcy_to_m2;
      
        // bool darcy_filter = kr_val*settings_->darcy_perm > perm_threshold;
      
        filter_functor<phase1> filter{pni_, &state_, false};
      
        auto [nrows, mapping] = pni_->generate_mapping(*settings_->solver.decomposition, filter);
      
        std::conditional_t<phase1, phase1_coef_functor, phase0_coef_functor> term{pni_, &state_, theta_adv, -999.999};
      
        auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);
      
        using namespace std::chrono;
      
        auto hash = std::hash<std::string_view>{}(
          std::string_view{reinterpret_cast<char*>(input.values.get()), nvalues*sizeof(HYPRE_Complex)});


        std::cout << fmt::format("ph{} | hash: {:x}", // kr: {:.2e}, dar_sw: {:.4f}, 
          int{phase1}, hash); /*kr_val, darcy_sw, */
      
      
        
      
        if (auto found = pressure_cache::cache().find(hash);
          found == pressure_cache::cache().end()) {
          dpl::strong_vector<net_tag, double> pressure(pni_->connected_count());
          auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);
      
          {
            // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);
      
            auto t0 = high_resolution_clock::now();
            dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_->solver.tolerance, settings_->solver.max_iterations);
            auto t1 = high_resolution_clock::now();
            std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings_->solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
            auto t2 = high_resolution_clock::now();
            std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);
      
            std::cout << fmt::format(" | matrix save: {} s, solve: {} s\n",
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
        else {
          std::cout << '\n';
          return found->second;  
        }
      
        
      };
      
      auto rel_calc_report = [=, this](double sw) {
        auto def = rel_calc(std::false_type{});
        auto inv = rel_calc(std::true_type{});
      
        kr_curves_.emplace_back(sw, def/absolute_rate, inv/absolute_rate);
      };
















      auto inv_idx_ = 0;

      for (; !queue.empty(); ++progress_idx_) {
        // if (++inv_idx_ < 1000)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{30});



        if (auto sw = eval_inv_volume()/total_pore_volume_; (sw - last_kr_sw) > 0.075) {
          last_kr_sw = sw;
          rel_calc_report(sw);
        }




        if (auto pc_point = eval_pc_point();
          std::abs(last_pc_point_.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point_.y() - pc_point.y()) > pc_max_step_) {
          last_pc_point_ = pc_point;
          add_to_pc_curve(pc_point);
        }

        auto [elem, local_idx, r_cap] = queue.front();

        queue.pop();

        if (elem == displ_elem::macro) {
          auto macro_idx = macro_idx_t(local_idx);
          auto net_idx = pni_->net(macro_idx);

          if (
            et_algo::get_header(dc_props.get_entry(*net_idx)) !=
            et_algo::get_header(outlet_entry))
            continue;

          state_.r_cap_global = std::max(r_cap, state_.r_cap_global);
          state_.config(net_idx) = phase_config{phase_config::phase0()};
          
          auto vol = volume(pn_, macro_idx);
          area_corner_mult -= vol/eq_tr::area(r_ins(pn_, macro_idx));
          inv_volume_coef0 += vol;

          {
            auto v = *net_idx;

            for (auto de : graph_.edges(v))
              if (!dc_props.is_null_entry(de) && !dc_props.is_tree_edge(de))
                dc_context_.non_tree_edge_remove(de);

            for (auto de : graph_.edges(v))
              if (!dc_props.is_null_entry(de))
                if (!dc_context_.tree_edge_split_and_reconnect(de))
                  if (auto hdr = et_algo::get_header(dc_props.get_entry(target(de, graph_))); hdr != et_algo::get_header(outlet_entry))
                    freeze_cluster(hdr);

            dc_props.set_entry(v, nullptr);
          }

          
        }
        else if (elem == displ_elem::throat) {
          auto [l, r] = pn_->throat_[adj][local_idx];

          et_traits::node_ptr outlet_hdr = et_algo::get_header(outlet_entry);

          et_traits::node_ptr l_et = dc_props.get_entry(*pni_->net(l));
          et_traits::node_ptr r_et = r == pn_->inlet() ? nullptr : dc_props.get_entry(*pni_->net(r));
          
          if (!(
            (l_et && et_algo::get_header(l_et) == outlet_hdr) ||
            (r_et && et_algo::get_header(r_et) == outlet_hdr)))
            continue;

          state_.r_cap_global = std::max(r_cap, state_.r_cap_global);
          state_.config(local_idx) = phase_config{phase_config::phase0()};

          
          auto vol = volume(pn_, local_idx);
          area_corner_mult -= vol/eq_tr::area(r_ins(pn_, local_idx));
          inv_volume_coef0 += vol;

          if (auto de = throat_to_de_[local_idx]; !dc_props.is_null_entry(de))
            if (dc_props.is_tree_edge(de))
              if (!dc_context_.tree_edge_split_and_reconnect(de)) {
                auto hdr = et_algo::get_header(dc_props.get_entry(target(de, graph_)));
                if (hdr == et_algo::get_header(outlet_entry))
                  hdr = et_algo::get_header(dc_props.get_entry(target(opposite(de, graph_), graph_)));

                freeze_cluster(hdr);
              }
              else
                dc_context_.non_tree_edge_remove(de);
        }

      }



      // for (std::size_t i{0}; i < pn_->throat_count(); ++i)
      //   if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inlet() == r) // inlet macro nodes TODO: voxels inlet needed
      //     if (pni_->connected(l)) {
      //       queue.insert(displ_elem::throat, i, eq_tr::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][i]));
      //       explored_throat[i] = true;
      //     }
    }

    void launch_primary(
      double absolute_rate,
      double theta,
      std::span<const dpl::vector2d> darcy_pc_inv) {

      pressure_cache::load();

      using eq_tr = hydraulic_properties::equilateral_triangle_properties;

      double darcy_r_cap_const;

      if (darcy_pc_inv.empty()) {
        auto min_r_cap_throat = std::numeric_limits<double>::max();
      
        for (std::size_t i{0}; i < pn_->throat_count(); ++i)
          if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && pni_->connected(l))
            min_r_cap_throat = std::min(min_r_cap_throat, pn_->throat_[attribs::r_ins][i]);
        min_r_cap_throat = 0.95*eq_tr::r_cap_piston_with_films(theta, min_r_cap_throat);
      
        darcy_r_cap_const = min_r_cap_throat;
      }
      else
        darcy_r_cap_const = 1/darcy_pc_inv.front().x();

      pc_max_step_ = 0.075/darcy_r_cap_const;

      dpl::strong_vector<net_tag, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue<true> queue;

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inlet() == r) // inlet macro nodes TODO: voxels inlet needed
          if (pni_->connected(l)) {
            queue.insert(displ_elem::throat, i, eq_tr::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][i]));
            explored_throat[i] = true;
          }

      net_idx_t outlet_idx{*pni_->connected_count()};
      dpl::graph::dc_properties dc_props{graph_};
      dpl::graph::et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);

      dpl::vector3d inv_volume_coefs{0};
      idx1d_t inv_darcy_count{0};

      state_.r_cap_global = queue.front().radius_cap;

      

      last_pc_point_ = {2, -1};
      double last_kr_sw = 1.0;
      kr_curves_.emplace_back(1, 1, 0);

      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim()).prod()*settings_->darcy_poro;

      auto eval_inv_volume = [&] {
        auto macro = inv_volume_coefs[0] + inv_volume_coefs[2]*state_.r_cap_global*state_.r_cap_global;
        if (!darcy_pc_inv.empty())
          macro += (1 - solve(darcy_pc_inv, 1/state_.r_cap_global, dpl::extrapolant::flat))*unit_darcy_pore_volume*inv_darcy_count;
        return macro;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{1 - eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };

      auto rel_calc = [=, this](auto phase1) {
        auto darcy_sw = settings_->darcy_kr0.empty()
          ? std::numeric_limits<double>::quiet_NaN()
          : solve(darcy_pc_inv, 1/state_.r_cap_global, dpl::extrapolant::flat);

        double kr_val = 0.0;

        if constexpr (phase1) {
          if (!settings_->darcy_kr1.empty())
            kr_val = std::max(0.0, solve(std::span{settings_->darcy_kr1}, darcy_sw, dpl::extrapolant::flat));
        }
        else {
          if (!settings_->darcy_kr0.empty())
            kr_val = std::max(0.0, solve(std::span{settings_->darcy_kr0}, darcy_sw, dpl::extrapolant::flat));
        }

        static constexpr auto kr_threshold = 1e-3;

        bool darcy_filter = kr_val > kr_threshold;

        filter_functor<phase1> filter{pni_, &state_, darcy_filter};
        auto [nrows, mapping] = pni_->generate_mapping(*settings_->solver.decomposition, filter);

        std::conditional_t<phase1, phase1_coef_functor, phase0_coef_functor> term{pni_, &state_, theta, kr_val*settings_->darcy_perm};

        auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);

        using namespace std::chrono;

        auto hash = std::hash<std::string_view>{}(
          std::string_view{reinterpret_cast<char*>(input.values.get()), nvalues*sizeof(HYPRE_Complex)});

        std::cout << fmt::format("ph{} | K: {:.2e} mD | hash: {:x}", // kr: {:.2e}, dar_sw: {:.4f}, 
          int{phase1}, kr_val*settings_->darcy_perm/presets::darcy_to_m2, hash); /*kr_val, darcy_sw, */


        auto found = pressure_cache::cache().find(hash);

        if (found == pressure_cache::cache().end()) {
          dpl::strong_vector<net_tag, double> pressure(pni_->connected_count());
          auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);

          {
            // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);

            auto t0 = high_resolution_clock::now();
            dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_->solver.tolerance, settings_->solver.max_iterations);
            auto t1 = high_resolution_clock::now();
            std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings_->solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
            auto t2 = high_resolution_clock::now();
            std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);

            std::cout << fmt::format(" | matrix save: {} s, solve: {} s\n",
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
      };

      auto rel_calc_report = [=, this](double sw) {
        auto def = rel_calc(std::false_type{});
        auto inv = rel_calc(std::true_type{});

        kr_curves_.emplace_back(sw, def/absolute_rate, inv/absolute_rate);
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
          add_to_pc_curve(pc_point);
        }

        auto [elem, local_idx, r_cap] = queue.front();

        queue.pop();

        if (elem == displ_elem::macro) {
          auto net_idx = pni_->net(macro_idx_t(local_idx));

          // {
          //   if (
          //     dpl::graph::et_algo::get_header(dc_props.get_entry(*net_idx)) !=
          //     dpl::graph::et_algo::get_header(outlet_entry))
          //     continue;
          //
          //   if (!props::has_films(theta))
          //     dc_context_.adjacent_edges_remove(*net_idx, dc_graph_); // TODO
          // }

          state_.config(net_idx) = phase_config::phase1_bulk_films(); // TODO

          state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

          inv_volume_coefs[0] += pn_->node_[attribs::volume][local_idx];
          inv_volume_coefs[2] += pn_->node_[attribs::volume][local_idx]*
            -eq_tr::area_films(theta)/eq_tr::area(pn_->node_[attribs::r_ins][local_idx]);


          for (auto ab : graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, graph_)}; b_net_idx != outlet_idx) // not outlet
              if (pni_->is_macro(b_net_idx)) { // macro (throat)
                if (auto throat_idx = de_to_throat_[ab]; !explored_throat[throat_idx]) {
                  queue.insert(displ_elem::throat, throat_idx, eq_tr::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]));
                  explored_throat[throat_idx] = true;
                }
              }
              else { // darcy
                if (!explored[b_net_idx]) {
                  auto b_voxel_idx = pni_->voxel(b_net_idx);
                  queue.insert(displ_elem::voxel, *b_voxel_idx, /*1/microprs_pc.back().y()*/ darcy_r_cap_const);
                  explored[b_net_idx] = true;
                }
              }
        }
        else if (elem == displ_elem::voxel) {
          auto net_idx = pni_->net(voxel_idx_t(local_idx));

          // {
          //   if (
          //     dpl::graph::et_algo::get_header(dc_props.get_entry(*net_idx)) !=
          //     dpl::graph::et_algo::get_header(outlet_entry))
          //     continue;
          //
          //   if (!props::has_films(theta))
          //     dc_context_.adjacent_edges_remove(*net_idx, dc_graph_); // TODO
          // }

          darcy_invaded_ = true;

          state_.config(net_idx) = phase_config::phase1_bulk_films(); // TODO

          state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

          ++inv_darcy_count;

          for (auto ab : graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, graph_)}; b_net_idx != outlet_idx) // not outlet
              if (!explored[b_net_idx])
                if (pni_->is_macro(b_net_idx)) { // macro
                  auto b_macro_idx = pni_->macro(b_net_idx);
                  queue.insert(displ_elem::macro, *b_macro_idx, eq_tr::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*b_macro_idx]));
                  explored[b_net_idx] = true;
                }
                else { // darcy
                  auto b_voxel_idx = pni_->voxel(b_net_idx);
                  queue.insert(displ_elem::voxel, *b_voxel_idx, /*1/microprs_pc.back().y()*/ darcy_r_cap_const);
                  explored[b_net_idx] = true;
                }
        }
        else { // throat
          auto [l, r] = pn_->throat_[attribs::adj][local_idx];

          // {
          //   using et_ptr = dpl::graph::et_traits::node_ptr;
          //
          //   et_ptr outlet_hdr = dpl::graph::et_algo::get_header(outlet_entry);
          //
          //   et_ptr l_entry = dc_props.get_entry(*pni_->net(l));
          //   et_ptr r_entry = r == pn_->inlet() ? nullptr : dc_props.get_entry(*pni_->net(r));
          //
          //   if (!(
          //       (l_entry && dpl::graph::et_algo::get_header(l_entry) == outlet_hdr) ||
          //       (r_entry && dpl::graph::et_algo::get_header(r_entry) == outlet_hdr)))
          //     continue;
          //
          //
          //   if (!props::has_films(theta))
          //     if (auto de = throat_to_de_[local_idx]; !dc_props.is_null_entry(de))
          //       if (dc_props.is_tree_edge(de))
          //         dc_context_.tree_edge_split_and_reconnect(de);
          //       else
          //         dc_context_.non_tree_edge_remove(de);
          // }

          state_.config(local_idx) = phase_config::phase1_bulk_films(); // TODO

          if (pn_->inner_node(r)/* && pni_->connected(l)*/) {
            state_.r_cap_global = std::min(r_cap, state_.r_cap_global);

            inv_volume_coefs[0] += 
              pn_->throat_[attribs::volume][local_idx];

            inv_volume_coefs[2] += 
              pn_->throat_[attribs::volume][local_idx]*
              -eq_tr::area_films(theta)/eq_tr::area(pn_->throat_[attribs::r_ins][local_idx]);


            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(displ_elem::macro, *r, eq_tr::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*r]));
              explored[r_net] = true;
            }
          }

          if (auto l_net = pni_->net(l); !explored[l_net]) {
            queue.insert(displ_elem::macro, *l, eq_tr::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*l]));
            explored[l_net] = true;
          }
        }
      }


      if (!darcy_pc_inv.empty()) {
        auto end_pc = 1/state_.r_cap_global;
        auto max_dacry_pc = darcy_pc_inv.back().x();

        add_to_pc_curve(eval_pc_point());

        { /* Capillary pressure */
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_dacry_pc) - std::log10(end_pc))/steps);

          for (auto i = 0; i < steps; ++i) {
            state_.r_cap_global /= step;
            add_to_pc_curve(eval_pc_point());
          }
        }

        ++progress_idx_;

        state_.r_cap_global = 1/end_pc;

        auto darcy_sw = solve(darcy_pc_inv, state_.r_cap_global, dpl::extrapolant::flat);

        if (!settings_->darcy_kr0.empty()) {
          auto steps = 6;

          for (auto i = 0; i < steps; ++i) {
            auto target_sw = darcy_sw*(1 - i*(1./steps));
            // std::cout << fmt::format("{:.5f}\n", target_sw);
            state_.r_cap_global = 1/solve(std::span{settings_->darcy_pc}, target_sw, dpl::extrapolant::flat);

            rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
            ++progress_idx_;
          }

          kr_curves_.emplace_back(0, 0, 1);
          ++progress_idx_;
        }
      }

      finished_ = true;

      {
        generate_graph();
        launch_secondary(absolute_rate, theta, theta);
      }
    }
  };
}