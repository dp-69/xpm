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

  class invasion_task {

    struct phase_state
    {
      dpl::strong_vector<net_tag, phase_config> macro_voxel;
      std::vector<phase_config> throat;

      double r_cap;
    };

    template <bool phase1>
    class filter_functor
    {
    public:
      const pore_network_image* pni;
      const phase_state* state;
      bool darcy_filter;

      bool operator()(macro_idx_t i) const {
        if constexpr (phase1)
          return state->macro_voxel[pni->net(i)].phase() == phase_config::phase1();
        else
          return true /*state->macro_voxel[pni->net(i)].phase() == phase_config::phase0()*/; // TODO
      }

      bool operator()(std::size_t i) const {
        if constexpr (phase1)
          return state->throat[i].phase() == phase_config::phase1();
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
      auto func(double area) const {
        using props = hydraulic_properties::equilateral_triangle_properties;
        return props::conductance_single_phase(area)*(1 - props::area_of_films(theta, state->r_cap)/area);
      }

    public:
      const pore_network_image* pni;
      const phase_state* state;
      double theta;
      double darcy_perm;
      

      double operator()(std::size_t i) const {
        using props = hydraulic_properties::equilateral_triangle_properties;
        return func(props::area(pni->pn().throat_[attribs::r_ins][i]));
      }

      double operator()(macro_idx_t i) const {
        using props = hydraulic_properties::equilateral_triangle_properties;
        return func(props::area(pni->pn().node_[attribs::r_ins][*i]));
      }

      double operator()(voxel_idx_t i) const {
        return darcy_perm; // TODO
      }
    };


    class phase0_coef_functor
    {
      auto cond_films() const {
        using props = hydraulic_properties::equilateral_triangle_properties;
        return props::conductance_films(theta, props::area_of_films(theta, state->r_cap));
      }

    public:
      const pore_network_image* pni;
      const phase_state* state;
      double theta;
      double darcy_perm;

      double operator()(std::size_t i) const {
        using props = hydraulic_properties::equilateral_triangle_properties;

        return
          state->throat[i].layout() == phase_config::bulk_films()
            ? cond_films()
            : props::conductance_single_phase(props::area(pni->pn().throat_[attribs::r_ins][i]));
      }

      double operator()(macro_idx_t i) const {
        using props = hydraulic_properties::equilateral_triangle_properties;

        return
          state->macro_voxel[pni->net(i)].layout() == phase_config::bulk_films()
            ? cond_films()
            : props::conductance_single_phase(props::area(pni->pn().node_[attribs::r_ins][*i]));
      }

      double operator()(voxel_idx_t i) const {
        return darcy_perm; // TODO
      }
    };






    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;

    dpl::graph::dc_graph dc_graph_;
    std::unordered_map<dpl::graph::dc_graph::edge_descriptor, std::size_t> de_to_throat_;
    std::unique_ptr<dpl::graph::dc_graph::edge_descriptor[]> throat_to_de_;

    dpl::graph::dc_context<dpl::graph::dc_properties> dc_context_;

    phase_state state_;

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
    explicit invasion_task(pore_network_image& pni) :
      pni_{&pni},
      pn_{&pni.pn()},
      img_{&pni.img()} {}

    auto& graph() {
      return dc_graph_;
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
      using clock = std::chrono::high_resolution_clock;
      using seconds = std::chrono::seconds;

      auto t0 = clock::now();
      std::cout << "graph...";

      std::tie(dc_graph_, de_to_throat_, throat_to_de_) = pni_->generate_dc_graph<true>();
      std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
        duration_cast<seconds>(clock::now() - t0).count(),
        dc_graph_.vertex_count(),
        dc_graph_.edge_count());

      auto t1 = clock::now();
      std::cout << "euler tour...";
      dc_context_.init_with_dfs(dc_graph_, dpl::graph::dc_properties{dc_graph_});
      std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t1).count());

      std::cout << "full decremental connectivity... (async)\n";

      pc_curve_.reserve(1000);
      kr_curves_.reserve(1000);
    }

    void launch_primary(
      const startup_settings& settings,
      double absolute_rate,
      double theta,
      std::span<const dpl::vector2d> darcy_pc_inv) {

      pressure_cache::load();

      using props = hydraulic_properties::equilateral_triangle_properties;

      state_.macro_voxel.resize(*pni_->connected_count());
      state_.throat.resize(pn_->throat_count());

      double darcy_r_cap_const;

      if (darcy_pc_inv.empty()) {
        auto min_r_cap_throat = std::numeric_limits<double>::max();
      
        for (std::size_t i{0}; i < pn_->throat_count(); ++i)
          if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && pni_->connected(l))
            min_r_cap_throat = std::min(min_r_cap_throat, pn_->throat_[attribs::r_ins][i]);
        min_r_cap_throat = 0.95*props::r_cap_piston_with_films(theta, min_r_cap_throat);
      
        darcy_r_cap_const = min_r_cap_throat;
      }
      else
        darcy_r_cap_const = 1/darcy_pc_inv.front().x();

      dpl::strong_vector<net_tag, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue<pressure_less> queue;

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inlet() == r) // inlet macro nodes TODO: voxels inlet needed
          if (pni_->connected(l)) {
            queue.insert(displ_elem::throat, i, props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][i]));
            explored_throat[i] = true;
          }

      net_idx_t outlet_idx{*pni_->connected_count()};

      dpl::graph::dc_properties dc_props{dc_graph_};
      dpl::graph::et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);

      dpl::vector3d inv_volume_coefs{0};
      idx1d_t inv_darcy_count{0};

      state_.r_cap = queue.front().radius_cap;

      auto total_pore_volume = pni_->total_pore_volume(settings.darcy_poro);

      dpl::vector2d last_pc_point{2, -1};
      double last_kr_sw = 1.0;
      kr_curves_.emplace_back(1, 1, 0);

      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim()).prod()*settings.darcy_poro;

      auto eval_inv_volume = [&] {
        return
          inv_volume_coefs[0] + inv_volume_coefs[2]*state_.r_cap*state_.r_cap
        + (1 - (darcy_pc_inv.empty() ? 0.0 : solve(darcy_pc_inv, 1/state_.r_cap, dpl::extrapolant::flat)))*unit_darcy_pore_volume*inv_darcy_count
        ;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{1 - eval_inv_volume()/total_pore_volume, 1/state_.r_cap}; };

      auto rel_calc = [=, this](auto phase1) {
        auto darcy_sw = settings.darcy_kr0.empty()
          ? std::numeric_limits<double>::quiet_NaN()
          : solve(darcy_pc_inv, 1/state_.r_cap, dpl::extrapolant::flat);

        double kr_val = 0.0;

        if constexpr (phase1) {
          if (!settings.darcy_kr1.empty())
            kr_val = std::max(0.0, solve(std::span{settings.darcy_kr1}, darcy_sw, dpl::extrapolant::flat));
        }
        else {
          if (!settings.darcy_kr0.empty())
            kr_val = std::max(0.0, solve(std::span{settings.darcy_kr0}, darcy_sw, dpl::extrapolant::flat));
        }

        static constexpr auto perm_threshold = 1e-6*presets::darcy_to_m2;

        bool darcy_filter = kr_val*settings.darcy_perm > perm_threshold;


        filter_functor<phase1> filter{pni_, &state_, darcy_filter};
        auto [nrows, mapping] = pni_->generate_mapping(*settings.solver.decomposition, filter);

        std::conditional_t<phase1, phase1_coef_functor, phase0_coef_functor> term{pni_, &state_, theta, kr_val*settings.darcy_perm};

        auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);

        using namespace std::chrono;

        auto hash = std::hash<std::string_view>{}(
          std::string_view{reinterpret_cast<char*>(input.values.get()), nvalues*sizeof(HYPRE_Complex)});

        std::cout << fmt::format("ph{} | K: {:.2e} mD | hash: {:x}\n", // kr: {:.2e}, dar_sw: {:.4f}, 
          int{phase1}, kr_val*settings.darcy_perm/presets::darcy_to_m2, hash); /*kr_val, darcy_sw, */


        auto found = pressure_cache::cache().find(hash);

        if (found == pressure_cache::cache().end()) {
          dpl::strong_vector<net_tag, double> pressure(*pni_->connected_count());
          auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);

          {
            // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings.solver.tolerance, settings.solver.max_iterations);

            auto t0 = high_resolution_clock::now();
            dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings.solver.tolerance, settings.solver.max_iterations);
            auto t1 = high_resolution_clock::now();
            std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings.solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
            auto t2 = high_resolution_clock::now();
            std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);

            std::cout << fmt::format("save matrix: {} s, solve matrix: {} s\n",
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

        if (auto sw = (1 - eval_inv_volume()/total_pore_volume); (last_kr_sw - sw) > 0.075) {
          last_kr_sw = sw;
          rel_calc_report(sw);
        }

        if (auto pc_point = eval_pc_point();
          std::abs(last_pc_point.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point.y() - pc_point.y()) > 0.075/darcy_r_cap_const) {
          last_pc_point = pc_point;
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

          state_.macro_voxel[net_idx] = phase_config::phase1_bulk_films(); // TODO

          state_.r_cap = std::min(r_cap, state_.r_cap);

          inv_volume_coefs[0] += pn_->node_[attribs::volume][local_idx];
          inv_volume_coefs[2] += pn_->node_[attribs::volume][local_idx]*
            -props::area_of_films(theta)/props::area(pn_->node_[attribs::r_ins][local_idx]);


          for (auto ab : dc_graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) // not outlet
              if (pni_->is_macro(b_net_idx)) { // macro (throat)
                if (auto throat_idx = de_to_throat_[ab]; !explored_throat[throat_idx]) {
                  queue.insert(displ_elem::throat, throat_idx, props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]));
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

          state_.macro_voxel[net_idx] = phase_config::phase1_bulk_films(); // TODO

          state_.r_cap = std::min(r_cap, state_.r_cap);

          ++inv_darcy_count;

          for (auto ab : dc_graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) // not outlet
              if (!explored[b_net_idx])
                if (pni_->is_macro(b_net_idx)) { // macro
                  auto b_macro_idx = pni_->macro(b_net_idx);
                  queue.insert(displ_elem::macro, *b_macro_idx, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*b_macro_idx]));
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

          state_.throat[local_idx] = phase_config::phase1_bulk_films(); // TODO

          if (pn_->inner_node(r)/* && pni_->connected(l)*/) {
            state_.r_cap = std::min(r_cap, state_.r_cap);

            inv_volume_coefs[0] += 
              pn_->throat_[attribs::volume][local_idx];

            inv_volume_coefs[2] += 
              pn_->throat_[attribs::volume][local_idx]*
              -props::area_of_films(theta)/props::area(pn_->throat_[attribs::r_ins][local_idx]);


            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(displ_elem::macro, *r, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*r]));
              explored[r_net] = true;
            }
          }

          if (auto l_net = pni_->net(l); !explored[l_net]) {
            queue.insert(displ_elem::macro, *l, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*l]));
            explored[l_net] = true;
          }
        }
      }


      if (!darcy_pc_inv.empty()) {
        auto end_pc = 1/state_.r_cap;
        auto max_dacry_pc = darcy_pc_inv.back().x();

        add_to_pc_curve(eval_pc_point());

        { /* Capillary pressure */
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_dacry_pc) - std::log10(end_pc))/steps);

          for (auto i = 0; i < steps; ++i) {
            state_.r_cap /= step;
            add_to_pc_curve(eval_pc_point());
          }
        }

        ++progress_idx_;

        state_.r_cap = 1/end_pc;

        auto darcy_sw = solve(darcy_pc_inv, state_.r_cap, dpl::extrapolant::flat);

        if (!settings.darcy_kr0.empty()) {
          auto steps = 6;

          for (auto i = 0; i < steps; ++i) {
            auto target_sw = darcy_sw*(1 - i*(1./steps));
            // std::cout << fmt::format("{:.5f}\n", target_sw);
            state_.r_cap = 1/solve(std::span{settings.darcy_pc}, target_sw, dpl::extrapolant::flat);

            rel_calc_report(1 - eval_inv_volume()/total_pore_volume);
            ++progress_idx_;
          }

          kr_curves_.emplace_back(0, 0, 1);
          ++progress_idx_;
        }
      }

      finished_ = true;
    }
  };
}