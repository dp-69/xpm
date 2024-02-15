#pragma once

#include "displ_queue.hpp"
#include "pore_network_image.hpp"

#include <dpl/graph/dc_context.hpp>

namespace xpm {
  

  class phase_state
  {
    dpl::strong_vector<net_t, phase_config> config_;
    std::vector<phase_config> config_t_;

    dpl::strong_vector<net_t, double> local_;
    std::vector<double> local_t_;

    static inline constexpr auto r_cap_mobile_ = std::numeric_limits<double>::min();

  public:
    double r_cap_global;

    void resize(net_t net_count, std::size_t throat_count) {
      config_.resize(net_count);
      local_.resize(net_count, r_cap_mobile_);

      config_t_.resize(throat_count);
      local_t_.resize(throat_count, r_cap_mobile_);
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

    void write_occupancy_image(
      const std::filesystem::path& path,
      const pore_network_image& pni,
      const double theta,
      const std::span<const dpl::vector2d> pc_inv) const {

      const auto& img = pni.img();
      const auto& dim = img.dim();
      const auto& pn = pni.pn();
      using type = std::uint16_t;
      static constexpr auto max = std::numeric_limits<type>::max();

      type darcy_value = pc_inv.empty() ? 0 : (1 - solve(pc_inv, 1/r_cap_global, dpl::extrapolant::flat))*max;

      using namespace presets;
      using namespace attrib;
      using eq_tr = hydraulic_properties::equilateral_triangle;

      std::vector<type> output(*img.size());

      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      voxel_t idx1d{0};

      for (k = 0; k < dim.z(); ++k)
        for (j = 0; j < dim.y(); ++j)
          for (i = 0; i < dim.x(); ++i, ++idx1d)
            if (img.dict.is_void(img.phase[idx1d])) {
              if (auto velem = img.velem[idx1d]; velem)
                if (macro_t macro{velem}; pni.connected(macro))
                  if (auto net = pni.net(macro); config(net).phase() == phase_config::phase1())
                    output[*idx1d] = (1 - eq_tr::area_corners(theta, r_cap(net))/eq_tr::area(r_ins(pn, macro)))*max;  // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            else if (img.dict.is_darcy(img.phase[idx1d]))
              output[*idx1d] = darcy_value;

      std::ofstream{path, std::ios::binary}
        .write(reinterpret_cast<char*>(output.data()), sizeof(type)**img.size());
    }
  };

  class invasion_task {
    template <bool phase1>
    struct filter_phase
    {
      const pore_network_image* pni;
      const phase_state* state;
      bool darcy_filter;

      bool operator()(macro_t i) const {
        if constexpr (phase1)
          return state->config(pni->net(i)).phase() == phase_config::phase1();
        else
          return true; // for theta < Pi/3 in eq. tr.
      }

      bool operator()(std::size_t i) const {
        if constexpr (phase1)
          return state->config(i).phase() == phase_config::phase1();
        else
          return true;
      }

      bool operator()(voxel_t) const {
        return darcy_filter;
      }
    };

    template <bool phase1, typename Perm/* = dpl::default_map::nan_t*/>
    struct coef_phase
    {
      const pore_network_image* pni;
      const phase_state* state;
      double theta;
      Perm darcy_perm;

      coef_phase(const pore_network_image* pni, const phase_state* state, double theta, Perm darcy_perm, std::bool_constant<phase1> = {})
        : pni(pni), state(state), theta(theta), darcy_perm(darcy_perm) {}

      auto operator()(macro_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;

        if constexpr (phase1) {
          auto area = eq_tr::area(attrib::r_ins(pni->pn(), i));
          return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(pni->net(i)))/area);
        }
        else
          return state->config(pni->net(i)).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(pni->net(i))))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni->pn(), i)));
      }

      auto operator()(std::size_t i) const {
        using eq_tr = hydraulic_properties::equilateral_triangle;
        if constexpr (phase1) {
          auto area = eq_tr::area(attrib::r_ins(pni->pn(), i));
          return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta, state->r_cap(i))/area);
        }
        else
          return state->config(i).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta, eq_tr::area_films(theta, state->r_cap(i)))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni->pn(), i)));
      }

      auto operator()(voxel_t i) const {
        return darcy_perm(i);
      }
    };

    using vertex_t = dpl::graph::dc_graph::vertex_t;
    using edge_t = dpl::graph::dc_graph::edge_t;
    using et_ptr = dpl::graph::et_traits::node_ptr;

    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;
    const runtime_settings* settings_;

    double total_pore_volume_;

    dpl::graph::dc_graph g_;
    std::unordered_map<edge_t, std::size_t> de_to_throat_;
    std::unique_ptr<edge_t[]> throat_to_de_;

    dpl::graph::dc_context context_;

    phase_state state_;
    std::vector<double> kr_recorded_;


    dpl::vector2d last_pc_point_;
    double pc_max_step_;

    struct curves
    {
      std::vector<dpl::vector2d> pc;
      std::vector<dpl::vector3d> kr;

      void add_pc(const dpl::vector2d& p) {
        pc.push_back(p);
      }
    } primary_,
      secondary_;

    std::filesystem::path image_dir_;

    idx1d_t progress_idx_ = 0;

    bool darcy_invaded_ = false;
    bool finished_ = false;

    

    

  public:
    explicit invasion_task(pore_network_image& pni, const runtime_settings& settings) :
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

    auto& primary() const {
      return primary_;
    }

    auto& secondary() const {
      return secondary_;
    }

    auto progress_idx() const {
      return progress_idx_;
    }

    bool finished() const {
      return finished_;
    }

    void init() {
      if (settings_->occupancy_images) {
        image_dir_ =
          std::filesystem::path(dpl::hypre::mpi::mpi_exec)
            .replace_filename("results")/settings_->image.path.stem()/"images";
          // std::filesystem::path(dpl::hypre::mpi::mpi_exec).replace_filename("image")/settings_->image.path.stem();

        remove_all(image_dir_);
        create_directories(image_dir_);
      }

      state_.resize(pni_->connected_count(), pn_->throat_count());

      primary_.pc.reserve(500);
      primary_.kr.reserve(500);
      secondary_.pc.reserve(500);
      secondary_.kr.reserve(500);

      total_pore_volume_ = pni_->total_pore_volume(
        [this](voxel_t i) { return settings_->darcy.poro(img_->phase[i]); });


      {
        using namespace std::chrono;
        std::cout << "graph...";
        auto t0 = high_resolution_clock::now();
        std::tie(g_, de_to_throat_, throat_to_de_) = pni_->generate_dc_graph();
        std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
          duration_cast<seconds>(high_resolution_clock::now() - t0).count(), g_.vertex_count(), g_.edge_count());
      }

      std::cout << "full decremental connectivity... (async)\n\n";
    }

    void generate_euler_tour() {
      using namespace std::chrono;
      using clock = high_resolution_clock;

      std::cout << "euler tour...";
      auto t0 = clock::now();

      context_.init_with_dfs(g_);
        
      for (std::size_t t{0}; t < pn_->throat_count(); ++t)
        if (auto [l, r] = attrib::adj(pn_, t); pni_->connected(l) && state_.config(t).phase() == phase_config::phase0())
          if (edge_t de = throat_to_de_[t]; !is_null_entry(de, g_) && !is_tree_edge(de, g_))
            context_.non_tree_edge_remove(de);

      for (std::size_t t{0}; t < pn_->throat_count(); ++t)
        if (auto [l, r] = attrib::adj(pn_, t); pni_->connected(l) && state_.config(t).phase() == phase_config::phase0())
          if (edge_t de = throat_to_de_[t]; !is_null_entry(de, g_))
            context_.tree_edge_split_and_reconnect(de);

      for (vertex_t v{0}; v < *pni_->connected_count(); ++v)
        if (net_t net{*v}; state_.config(net).phase() == phase_config::phase0())
          context_.adjacent_edges_remove(v, g_);

      std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
    }

    double calc_relative(const runtime_settings::input_curves& curve, std::span<const dpl::vector2d> pc_to_sw, double theta, auto phase1) {
      using dpl::extrapolant::flat;

      auto kr = pc_to_sw.empty() ? 0.0 : solve(curve.kr[phase1], solve(pc_to_sw, 1/state_.r_cap_global, flat), flat);

      // static constexpr auto kr_threshold = 1e-3;
      // bool darcy_filter = kr > kr_threshold;

      filter_phase<phase1> filter{pni_, &state_, kr > 1e-3/*darcy_filter*/};
      coef_phase term{pni_, &state_, theta, [kr, this](voxel_t i) { return kr*settings_->darcy.perm(img_->phase[i]); }, phase1};

      auto [nrows, mapping] = pni_->generate_mapping(*settings_->solver.decomposition, filter);
      auto [nvalues, input] = pni_->generate_pressure_input(nrows, mapping.forward, term, filter);

      auto hash = pressure_cache::hash(nvalues, input);

      fmt::print("ph{} | hash: {:x}", int{phase1}, hash);

      auto found = pressure_cache::cache().find(hash);

      if (!settings_->solver.cache.use || found == pressure_cache::cache().end()) {
        dpl::strong_vector<net_t, double> pressure(pni_->connected_count());
        auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);

        {
          using namespace std::chrono;

          // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);

          auto t0 = high_resolution_clock::now();
          dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_->solver.tolerance, settings_->solver.max_iterations);
          auto t1 = high_resolution_clock::now();
          std::system(fmt::format("mpiexec -np {} \"{}\" -s", settings_->solver.decomposition->prod(), dpl::hypre::mpi::mpi_exec).c_str()); // NOLINT(concurrency-mt-unsafe)
          auto t2 = high_resolution_clock::now();
          std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::mpi::load_values(nrows);

          std::cout << fmt::format(" | store: {} s, solve: {} s\n",
            duration_cast<seconds>(t1 - t0).count(), duration_cast<seconds>(t2 - t1).count());
        }
        
        for (HYPRE_BigInt i = 0; i < nrows; ++i)
          pressure[mapping.backward[i]] = decomposed_pressure[i];

        auto inlet = pni_->flow_rates(pressure, term, filter).second;

        if (settings_->solver.cache.save) {
          pressure_cache::cache()[hash] = inlet;
          pressure_cache::save();
        }

        return inlet;
      }

      std::cout << '\n';
      return found->second;
    }

    void launch_secondary(double absolute_rate, double theta) {
      using eq_tr = hydraulic_properties::equilateral_triangle;
      using namespace attrib;

      auto cell_volume = (pn_->physical_size/img_->dim()).prod();
      auto pc_to_sw = settings_->secondary.calc_pc_inv();
      auto darcy_r_cap = pc_to_sw.empty() ? std::numeric_limits<double>::quiet_NaN() : 1/settings_->secondary.pc.back().y();

      dpl::strong_vector<net_t, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue<false> queue;

      auto area_corner_mult = 0.0;
      auto inv_volume_coef0 = 0.0;


      auto bulk_films_filter = [this](auto i) { return state_.config(i).layout() == phase_config::bulk_films(); };

      for (macro_t i{0}; i < pn_->node_count(); ++i)
        if (pni_->connected(i)) {
          if (bulk_films_filter(pni_->net(i))) {
            area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));
            queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
          }
          else
            inv_volume_coef0 += volume(pn_, i);
        }
      
      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = adj(pn_, i); pni_->connected(l)) {
          if (bulk_films_filter(i)) {
            area_corner_mult += volume(pn_, i)/eq_tr::area(r_ins(pn_, i));

            if (pn_->inlet() == r) {
              queue.insert(i, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, i)));
              explored_throat[i] = true;
            }
            else
              queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
          }
          else
            inv_volume_coef0 += volume(pn_, i);
        }

      auto inv_total_porosity = 0.0;

      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_t idx1d{0};

        for (k = 0; k < img_->dim().z(); ++k)
          for (j = 0; j < img_->dim().y(); ++j)
            for (i = 0; i < img_->dim().x(); ++i, ++idx1d)
              if (pni_->connected(idx1d)) { // TODO: Early primary termination
                inv_total_porosity += settings_->darcy.poro(img_->phase[idx1d]);
                queue.insert(idx1d, darcy_r_cap);
              }
      }


      auto eval_inv_volume = [&] {
        auto macro = area_corner_mult*eq_tr::area_corners(theta, state_.r_cap_global) + inv_volume_coef0;
        if (!pc_to_sw.empty())
          macro += solve(pc_to_sw, 1/state_.r_cap_global, dpl::extrapolant::flat)*cell_volume*inv_total_porosity;
        return macro;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };

      using namespace dpl::graph;

      net_t outlet_idx = pni_->connected_count();
      et_ptr outlet_entry = get_entry(vertex_t(*outlet_idx), g_); 

      auto freeze_cluster = [&](et_ptr hdr) {
        auto area_corners = eq_tr::area_corners(theta, state_.r_cap_global);

        for (et_ptr et : range<et_traits>(hdr))
          if (is_loop_edge(et, g_)) {
            vertex_t v = get_vertex(et, g_);

            set_entry(v, nullptr, g_);

            net_t v_net{*v}; 

            if (bulk_films_filter(v_net)) {  // NOLINT(clang-diagnostic-float-equal)
              state_.local(v_net) = state_.r_cap_global;

              if (pni_->is_macro(v_net)) {
                auto macro_idx = pni_->macro(v_net);

                auto mult = volume(pn_, macro_idx)/eq_tr::area(r_ins(pn_, macro_idx));
                area_corner_mult -= mult;
                inv_volume_coef0 += area_corners*mult;  
              }
              else {
                auto voxel = pni_->voxel(v_net);
                auto poro = settings_->darcy.poro(img_->phase[voxel]);
                inv_total_porosity -= poro;
                inv_volume_coef0 += solve(pc_to_sw, 1/state_.r_cap_global, dpl::extrapolant::flat)*cell_volume*poro;  
              }
            }

            if (pni_->is_macro(v_net))
              for (edge_t vu : g_.edges(v))
                if (pni_->is_macro(net_t{*target(vu, g_)})) {
                  auto t_idx = de_to_throat_[vu];

                  if (bulk_films_filter(t_idx) && state_.mobile(t_idx)) {
                    set_null_entry(vu, g_);
                    set_null_entry(opposite(vu, g_), g_);  

                    state_.local(t_idx) = state_.r_cap_global;

                    auto mult = volume(pn_, t_idx)/eq_tr::area(r_ins(pn_, t_idx));
                    area_corner_mult -= mult;
                    inv_volume_coef0 += area_corners*mult;
                  }
                } 
                else
                  set_null_entry(vu, g_);
            else
              for (edge_t vu : g_.edges(v))
                set_null_entry(vu, g_);
          }
      };



      double last_kr_sw = 0.0;
      secondary_.kr.emplace_back(0, 0, 1);
      // secondary_.add_pc_point(last_pc_point_); // TODO




      auto rel_calc_report = [&](double sw) {
        auto phase1_connected = [&] {
          for (auto [l, r] : pn_->throat_.span(adj))
            if (r == pn_->inlet() && pni_->connected(l)) // macro-inlet
              if (get_entry(vertex_t(*pni_->net(l)), g_))
                return true;

          {
            idx3d_t ijk;
            auto& [i, j, k] = ijk;
          
            for (k = 0; k < img_->dim().z(); ++k)
              for (j = 0; j < img_->dim().y(); ++j)
                if (voxel_t idx1d{img_->idx_map(0, j, k)}; pni_->connected(idx1d)/*img_->phase[idx1d] == presets::microporous*/) // darcy-inlet
                  if (get_entry(vertex_t(*pni_->net(idx1d)), g_))
                    return true;
          }

          return false;
        };


        
        auto def = calc_relative(settings_->secondary, pc_to_sw, theta, std::false_type{});
        auto inv = phase1_connected() ? calc_relative(settings_->secondary, pc_to_sw, theta, std::true_type{}) : 0;
      
        secondary_.kr.emplace_back(sw, def/absolute_rate, inv/absolute_rate);
      };



      auto write_occupancy_image = [&](double sw) {
        if (settings_->occupancy_images)
          state_.write_occupancy_image(
            image_dir_/fmt::format("secondary-{:.4f}.raw", sw), *pni_, theta, pc_to_sw);
      };

      if (pc_to_sw.empty()) {
        last_pc_point_ = eval_pc_point();
        secondary_.add_pc(last_pc_point_); // TODO: simplify with else
      }
      else {
        auto max_pc = primary_.pc.back().y();

        state_.r_cap_global = 1/max_pc;
        secondary_.add_pc(eval_pc_point());
        
        {
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_pc) - std::log10(1/darcy_r_cap))/steps);

          for (auto i = 0; i < steps; ++i) {
            state_.r_cap_global *= step;
            last_pc_point_ = eval_pc_point();
            secondary_.add_pc(last_pc_point_);
            if (i%2 != 0)
              write_occupancy_image(last_pc_point_.x());
          }
        }

        ++progress_idx_;

        if (!settings_->secondary.kr[0].empty()) {
          auto steps = 6;

          for (auto i = 1; i < steps; ++i) {
            state_.r_cap_global = 1/solve(settings_->secondary.pc, 1.*i/steps, dpl::extrapolant::flat);
            auto sw = eval_inv_volume()/total_pore_volume_;
            last_kr_sw = sw;
            rel_calc_report(sw);
            ++progress_idx_;
          }
        }
      }


      

      auto progress_percolation = [&](auto idx, double r_cap) {
        if (r_cap > state_.r_cap_global) {
          if (auto pc_point = eval_pc_point();
            std::abs(last_pc_point_.x() - pc_point.x()) > /*0.025*/ settings_->report.sw_pc ||
            std::abs(last_pc_point_.y() - pc_point.y()) > pc_max_step_) {
            last_pc_point_ = pc_point;
            secondary_.add_pc(pc_point);
            write_occupancy_image(pc_point.x());
          }

          if (auto sw = eval_inv_volume()/total_pore_volume_; sw - last_kr_sw > /*0.025*/ settings_->report.sw_kr) {
            last_kr_sw = sw;
            rel_calc_report(sw);
          }

          state_.r_cap_global = r_cap;
        }

        state_.config(idx) = phase_config{phase_config::phase0()};
      };


      auto local_prog_idx = 0;

      for (; !queue.empty(); ++progress_idx_) {
        // {
        //   using namespace std::this_thread;
        //   using namespace std::chrono;
        //
        //   if (local_prog_idx < 1000) {
        //     if (local_prog_idx%20 == 0)
        //       sleep_for(milliseconds{40});
        //   }
        //   else if (local_prog_idx < 4000) {
        //     if (local_prog_idx%40 == 0)
        //       sleep_for(milliseconds{50});
        //   }
        //
        //   ++local_prog_idx;
        // }


        auto [elem, local, r_cap] = queue.front();

        queue.pop();

        if (elem == displ_elem::macro) {
          macro_t macro(local);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net = pni_->net(macro);

          if (state_.config(net).phase() == phase_config::phase0())
            continue;

          if (get_entry(vertex_t(*net), g_)) {
            progress_percolation(net, r_cap);

            area_corner_mult -= volume(pn_, macro)/eq_tr::area(r_ins(pn_, macro));
            inv_volume_coef0 += volume(pn_, macro);

            {
              vertex_t v(*net);

              for (edge_t vu : g_.edges(v))
                if (!is_null_entry(vu, g_) && !is_tree_edge(vu, g_))
                  context_.non_tree_edge_remove(vu);

              for (edge_t vu : g_.edges(v)) {
                if (net_t u_net_idx{*target(vu, g_)}; u_net_idx == outlet_idx || pni_->is_macro(u_net_idx)) { // macro-macro
                  if (auto t_idx = de_to_throat_[vu]; !explored_throat[t_idx]) {
                    queue.insert(t_idx, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, t_idx)));
                    explored_throat[t_idx] = true;
                  }
                }
                // else { // macro-darcy
                // }

                if (!is_null_entry(vu, g_))
                  if (!context_.tree_edge_split_and_reconnect(vu))
                    if (et_ptr hdr = et_algo::get_header(get_entry(target(vu, g_), g_));
                      hdr != et_algo::get_header(outlet_entry))
                      freeze_cluster(hdr);
              }

              set_entry(v, nullptr, g_);
            }
          }
        }
        else if (elem == displ_elem::throat) {
          if (state_.config(local).phase() == phase_config::phase0())
            continue;

          auto [l, r] = adj(pn_, local);

          auto l_net = pni_->net(l);

          et_ptr l_entry = get_entry(vertex_t(*l_net), g_);
          et_ptr r_entry =
            r == pn_->inlet() ? nullptr :
            r == pn_->outlet() ? outlet_entry :
            get_entry(vertex_t(*pni_->net(r)), g_);

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

            progress_percolation(local, r_cap);

            area_corner_mult -= volume(pn_, local)/eq_tr::area(r_ins(pn_, local));
            inv_volume_coef0 += volume(pn_, local);

            if (edge_t de = throat_to_de_[local]; !is_null_entry(de, g_))
              if (is_tree_edge(de, g_)) {
                if (!context_.tree_edge_split_and_reconnect(de)) {
                  et_ptr hdr = et_algo::get_header(get_entry(target(de, g_), g_));
                  if (hdr == et_algo::get_header(outlet_entry))
                    hdr = et_algo::get_header(get_entry(target(opposite(de, g_), g_), g_));

                  freeze_cluster(hdr);
                }
              }
              else
                context_.non_tree_edge_remove(de);
          }
        }
        else { // voxel
          auto net = pni_->net(voxel_t(local));  // NOLINT(cppcoreguidelines-narrowing-conversions)
          
          // if (traits.get_entry(vertex_t(*net_idx))) 
          {
            progress_percolation(net, r_cap);

            vertex_t v(*net);
              
            for (edge_t vu : g_.edges(v))
              if (!is_null_entry(vu, g_) && !is_tree_edge(vu, g_))
                context_.non_tree_edge_remove(vu);
              
            for (edge_t vu : g_.edges(v))
              if (!is_null_entry(vu, g_) && !context_.tree_edge_split_and_reconnect(vu))
                if (et_ptr hdr = et_algo::get_header(get_entry(target(vu, g_), g_)); hdr != et_algo::get_header(outlet_entry))
                  freeze_cluster(hdr);
            
            set_entry(v, nullptr, g_);
          }
        }
      }


      secondary_.add_pc(eval_pc_point());
      rel_calc_report(eval_inv_volume()/total_pore_volume_);
      ++progress_idx_;
    }

    void launch_primary(double absolute_rate, double theta, const std::vector<dpl::vector2d>& pc_to_sw) {
      return;

      pressure_cache::load();

      // using namespace std::chrono;
      // auto t0 = high_resolution_clock::now();

      using eq_tr = hydraulic_properties::equilateral_triangle;
      using namespace attrib;
      
      auto darcy_r_cap = std::numeric_limits<double>::quiet_NaN();

      if (pc_to_sw.empty()) {
        auto min_r_cap_throat = std::numeric_limits<double>::max();

        for (std::size_t i{0}; i < pn_->throat_count(); ++i)
          if (auto [l, r] = adj(pn_, i); pn_->inner_node(r) && pni_->connected(l))
            min_r_cap_throat = std::min(min_r_cap_throat, r_ins(pn_, i));

        min_r_cap_throat = 0.95*eq_tr::r_cap_piston_with_films(theta, min_r_cap_throat);pc_max_step_ = 0.075/min_r_cap_throat;
      }
      else {
        darcy_r_cap = 1/pc_to_sw.front().x();
        pc_max_step_ = 0.075/darcy_r_cap;
      }

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
      auto inv_total_porosity = 0.0;

      state_.r_cap_global = queue.front().radius_cap;

      last_pc_point_ = {2, -1};
      double last_kr_sw = 1.0;
      primary_.kr.emplace_back(1, 1, 0);

      auto cell_volume = (pn_->physical_size/img_->dim()).prod();

      auto eval_inv_volume = [&] {
        auto macro = inv_volume_coefs[0] + inv_volume_coefs[2]*state_.r_cap_global*state_.r_cap_global;
        if (!pc_to_sw.empty())
          macro += (1 - solve(pc_to_sw, 1/state_.r_cap_global, dpl::extrapolant::flat))*cell_volume*inv_total_porosity;
        return macro;
      };

      auto eval_pc_point = [&] { return dpl::vector2d{1 - eval_inv_volume()/total_pore_volume_, 1/state_.r_cap_global}; };

      auto rel_calc_report = [this, /*t0, */theta, absolute_rate, &pc_to_sw](double sw) {
        // fmt::print("{}ms\n", duration_cast<milliseconds>(high_resolution_clock::now() - t0).count());
        primary_.kr.emplace_back(sw,
          calc_relative(settings_->primary, pc_to_sw, theta, std::false_type{})/absolute_rate,
          calc_relative(settings_->primary, pc_to_sw, theta, std::true_type{})/absolute_rate); // TODO: calculate when breakthrough
      };

      auto write_occupancy_image = [&](double sw) {
        if (settings_->occupancy_images)
          state_.write_occupancy_image(
            image_dir_/fmt::format("primary-{:.4f}.raw", sw), *pni_, theta, pc_to_sw);
      };

      auto progress_percolation = [&](auto idx, double r_cap) {
        if (r_cap < state_.r_cap_global) {
          if (auto pc_point = eval_pc_point();
            std::abs(last_pc_point_.x() - pc_point.x()) > settings_->report.sw_pc ||
            std::abs(last_pc_point_.y() - pc_point.y()) > pc_max_step_) {
            last_pc_point_ = pc_point;
            primary_.add_pc(pc_point);
            write_occupancy_image(pc_point.x());

            
          }

          if (auto sw = 1 - eval_inv_volume()/total_pore_volume_; last_kr_sw - sw > settings_->report.sw_kr) {
            last_kr_sw = sw;
            rel_calc_report(sw);
          }

          state_.r_cap_global = r_cap;
        }

        state_.config(idx) = phase_config::phase1_bulk_films();
      };

      


      progress_idx_ = 0;
      // getchar();
      for (/*progress_idx_ = 0*/; !queue.empty(); ++progress_idx_) {
        // {
        //   using namespace std::this_thread;
        //   using namespace std::chrono;
        //
        //   if (progress_idx_ < 35)
        //     sleep_for(milliseconds{150});
        //   else if (progress_idx_ < 100)
        //     sleep_for(milliseconds{50});
        //   else if (progress_idx_ < 300 /*1000*/)
        //     sleep_for(milliseconds{25});
        // }

        auto [elem, local, r_cap] = queue.front();

        if (1/r_cap > settings_->max_pc) {
          last_pc_point_ = eval_pc_point();
          primary_.add_pc(last_pc_point_);
          write_occupancy_image(last_pc_point_.x());
          rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
          break;
        }

        queue.pop();

        if (elem == displ_elem::macro) {
          auto macro = macro_t(local);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net = pni_->net(macro);

          progress_percolation(net, r_cap);

          inv_volume_coefs[0] += volume(pn_, macro);
          inv_volume_coefs[2] -= volume(pn_, macro)*eq_tr::area_films(theta)/eq_tr::area(r_ins(pn_, macro));
          
          for (edge_t vu : g_.edges(vertex_t(*net)))
            if (net_t u_net_idx{*target(vu, g_)}; u_net_idx == outlet_idx || pni_->is_macro(u_net_idx)) { // macro-macro
              if (auto t_idx = de_to_throat_[vu]; !explored_throat[t_idx]) {
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
          auto voxel = voxel_t(local);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net = pni_->net(voxel);

          darcy_invaded_ = true;
          inv_total_porosity += settings_->darcy.poro(img_->phase[voxel]);

          progress_percolation(net, r_cap); // TODO: phase_config::phase1_bulk_films(); of a voxel

          for (auto vu : g_.edges(vertex_t(*net))) {
            if (net_t u_net{*target(vu, g_)}; u_net != outlet_idx && !explored[u_net]) {
              if (pni_->is_macro(u_net)) { // macro
                auto u_macro = pni_->macro(u_net);
                queue.insert(u_macro, eq_tr::r_cap_piston_with_films(theta, r_ins(pn_, u_macro)));
                explored[u_net] = true;
              }
              else { // darcy
                queue.insert(pni_->voxel(u_net), darcy_r_cap);
                explored[u_net] = true;
              }
            }
          }
        }
        else { // throat
          auto [l, r] = adj(pn_, local);

          progress_percolation(local, r_cap);

          inv_volume_coefs[0] += volume(pn_, local);
          inv_volume_coefs[2] -= volume(pn_, local)*eq_tr::area_films(theta)/eq_tr::area(r_ins(pn_, local));

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


      if (!pc_to_sw.empty()) {
        auto end_pc = 1/state_.r_cap_global;
        auto max_darcy_pc = pc_to_sw.back().x();

        if (auto pc_point = eval_pc_point(); primary_.pc.back().x() != pc_point.x()) {  // NOLINT(clang-diagnostic-float-equal)
          primary_.add_pc(pc_point);
          write_occupancy_image(pc_point.x());
        }

        { /* Capillary pressure */
          auto steps = 10;

          auto step = std::pow(10, (std::log10(max_darcy_pc) - std::log10(end_pc))/steps);

          // auto i = 0;
          // for (auto next_r_cap = state_.r_cap_global/step; 1/next_r_cap < settings_->max_pc; next_r_cap /= step) {  // NOLINT(cert-flp30-c)
          //   state_.r_cap_global /= step;
          //   primary_.add_pc(eval_pc_point());
          //
          //   if (i++%2 != 0)
          //     write_occupancy_image(primary_.pc.back().x());
          // }

          for (auto i = 0; i < steps; ++i) {
            auto next_r_cap = state_.r_cap_global/step;
            if (1/next_r_cap > settings_->max_pc) {
              state_.r_cap_global = 1/settings_->max_pc;
              primary_.add_pc(eval_pc_point());
            
              if (i%2 != 0)
                write_occupancy_image(primary_.pc.back().x());

              break;
            }

            state_.r_cap_global = next_r_cap;
            primary_.add_pc(eval_pc_point());
          
            if (i%2 != 0)
              write_occupancy_image(primary_.pc.back().x());
          }

          ++progress_idx_;

          state_.r_cap_global = 1/end_pc;
        }

        

        // if (!settings_->primary.kr[0].empty()) 
        {
          auto steps = 6;

          // auto i = 0;
          // while (true)
          //   if (
          //     auto next_r_cap = 1/solve(settings_->primary.pc, 1 - 1.*i++/steps, dpl::extrapolant::flat);
          //     1/next_r_cap < settings_->max_pc) {
          //
          //     rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
          //     ++progress_idx_;
          //   }
          //   else
          //     break;

          for (auto i = 0; i < steps; ++i) {
            auto next_r_cap = 1/solve(settings_->primary.pc, 1 - 1.*i/steps, dpl::extrapolant::flat);
            if (1/next_r_cap > settings_->max_pc)
              break;

            state_.r_cap_global = next_r_cap;
            rel_calc_report(1 - eval_inv_volume()/total_pore_volume_);
            ++progress_idx_;
          }
        }
      }

      primary_.kr.emplace_back(0, 0, 1);
      ++progress_idx_;

      // getchar();

      finished_ = true;

      {
        std::cout << '\n';
        generate_euler_tour();
        launch_secondary(absolute_rate, theta);
      }
    }
  };
}