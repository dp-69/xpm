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

#include "displ_queue.hpp"
#include "pore_network_image.hpp"

#include <dpl/graph/dc_context.hpp>

namespace xpm {

  class phase_state
  {
    dpl::so_uptr<net_t, phase_config> config_;
    std::vector<phase_config> config_t_;

    dpl::so_uptr<net_t, double> local_;
    std::vector<double> local_t_;

    static constexpr auto r_cap_mobile_ = std::numeric_limits<double>::min();

  public:
    double r_cap_global;

    auto pc_global() const {
      return 1/r_cap_global;
    }

    void resize(net_t net_count, std::size_t throat_count) {
      config_.resize(net_count);
      local_.assign(net_count, r_cap_mobile_);

      config_t_.resize(throat_count);
      local_t_.resize(throat_count, r_cap_mobile_);
    }

    auto& config(net_t i) const { return config_[i]; }
    auto& config(net_t i) { return config_[i]; }

    auto& config(throat_t i) const { return config_t_[i]; }
    auto& config(throat_t i) { return config_t_[i]; }

    auto& local(net_t i) const { return local_[i]; }
    auto& local(net_t i) { return local_[i]; }

    auto& local(throat_t i) const { return local_t_[i]; }
    auto& local(throat_t i) { return local_t_[i]; }

    bool mobile(auto i) const {
      return this->local(i) == r_cap_mobile_;  // NOLINT(clang-diagnostic-float-equal)
    }

    double r_cap(auto i) const {
      auto l = this->local(i);
      return l == r_cap_mobile_ ? r_cap_global : l;  // NOLINT(clang-diagnostic-float-equal)
    }

    template<int cycle>
    void write_occupancy_image(
      const std::filesystem::path& path,
      const pore_network_image& pni,
      double theta,
      dpl::so_span<voxel_ns::phase_t, const darcy_info> darcy) const
    {
      const auto& img = pni.img();
      const auto& dim = img.dim();

      using type = std::uint16_t;
      static constexpr auto max{std::numeric_limits<type>::max()};

      using namespace attrib;
      using eq_tr = hydraulic_properties::equilateral_triangle;
      using voxel_ns::phase_t;

      dpl::so_uptr<phase_t, type> phase_output(darcy.size());
      for (phase_t i{0}; i < darcy.size(); ++i)
        phase_output[i] = (1 - darcy[i].pc_to_sw[cycle](pc_global()))*max;  // NOLINT(cppcoreguidelines-narrowing-conversions)

      std::vector<type> output(*img.size());

      idx3d_t ijk;
      auto& [i, j, k] = ijk;
      voxel_t idx1d{0};

      for (k = 0; k < dim.z; ++k)
        for (j = 0; j < dim.y; ++j)
          for (i = 0; i < dim.x; ++i, ++idx1d)
            if (img.is_void(idx1d)) {
              if (auto velem = img.velem[idx1d]; velem)
                if (macro_t macro{*velem}; pni.connected(macro))
                  if (auto net = pni.net(macro); config(net).phase() == phase_config::phase1())
                    output[*idx1d] = (1 - eq_tr::area_corners(theta, r_cap(net))/eq_tr::area(r_ins(pni.pn(), macro)))*max;  // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            else if (img.is_darcy(idx1d))
              output[*idx1d] = phase_output[img.phase[idx1d]];

      std::ofstream{path, std::ios::binary}
        .write(reinterpret_cast<char*>(output.data()), sizeof(type)**img.size());  // NOLINT(cppcoreguidelines-narrowing-conversions)
    }
  };

  class invasion_task {
    using phase_t = voxel_ns::phase_t;

    template <bool phase1>
    struct filter_t
    {
      const pore_network_image* pni;
      // const image_data* img;

      
      // const dpl::so_uptr<phase_t, darcy_info>* darcy_list;
      const dpl::so_uptr<voxel_t, phase_t>* phase_map;
      const dpl::so_uptr<phase_t, double>* kr_map;

      const phase_state* state;
      
      // bool darcy_filter;

      static constexpr auto kr_threshold = 1e-3;

      bool operator()(macro_t i) const {
        if constexpr (phase1)
          return state->config(pni->net(i)).phase() == phase_config::phase1();
        else
          return true; // for theta < Pi/3 in eq. tr.
      }

      bool operator()(throat_t i) const {
        if constexpr (phase1)
          return state->config(i).phase() == phase_config::phase1();
        else
          return true;
      }

      bool operator()(voxel_t i) const {
        return (*kr_map)[(*phase_map)[i]] > kr_threshold; 

        //
        // (*darcy_list)[(*phase_map)[i]].kr
        //
        // pni->img().phase[i]

        // return darcy_filter;
      }
    };

    template <bool phase1, typename Perm>
    class term_multi_t
    {
      using eq_tr = hydraulic_properties::equilateral_triangle;

      const pore_network_image* pni_;
      const phase_state* state_;
      double theta_;
      Perm darcy_perm_;

    public:
      explicit term_multi_t(std::bool_constant<phase1>,
        const pore_network_image* pni, const phase_state* state, double theta, Perm darcy_perm)
        : pni_(pni), state_(state), theta_(theta), darcy_perm_(darcy_perm) {}

      auto operator()(macro_t i) const {
        if constexpr (phase1) {
          auto area = eq_tr::area(attrib::r_ins(pni_->pn(), i));
          return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta_, state_->r_cap(pni_->net(i)))/area);
        }
        else
          return state_->config(pni_->net(i)).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta_, eq_tr::area_corners(theta_, state_->r_cap(pni_->net(i))))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni_->pn(), i)));
      }

      auto operator()(throat_t i) const {
        if constexpr (phase1) {
          auto area = eq_tr::area(attrib::r_ins(pni_->pn(), i));
          return eq_tr::conductance_single(area)*(1 - eq_tr::area_corners(theta_, state_->r_cap(i))/area);
        }
        else
          return state_->config(i).layout() == phase_config::bulk_films()
            ? eq_tr::conductance_films(theta_, eq_tr::area_corners(theta_, state_->r_cap(i)))
            : eq_tr::conductance_single(eq_tr::area(attrib::r_ins(pni_->pn(), i)));
      }

      auto operator()(voxel_t i) const { // TODO: BUG: Use frozen r_cap when trapped!??? If it is trapped, it is not flowing?
        return darcy_perm_(i);
      }
    };

    using vertex_t = dpl::graph::dc_graph::vertex_t;
    using edge_t = dpl::graph::dc_graph::edge_t;
    using et_ptr = dpl::graph::et_traits::node_ptr;

    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;
    const runtime_settings* cfg_;

    double total_pore_volume_;

    dpl::graph::dc_graph g_;
    std::unordered_map<edge_t, throat_t> de_to_throat_;
    std::unique_ptr<edge_t[]> throat_to_de_;

    dpl::graph::dc_context context_;

    phase_state state_;
    std::vector<double> kr_recorded_;


    dpl::vector2d last_pc_point_;
    double pc_max_log_step_;

    struct curves
    {
      std::vector<dpl::vector2d> pc;
      std::vector<dpl::vector3d> kr;
    } primary_,
      secondary_;

    std::filesystem::path image_dir_;

    idx1d_t progress_idx_ = 0;

    bool darcy_invaded_ = false;
    bool primary_finished_ = false;

    

    

  public:
    explicit invasion_task(pore_network_image& pni, const runtime_settings& cfg) :
      pni_{&pni},
      pn_{&pni.pn()},
      img_{&pni.img()},
      cfg_{&cfg} {}

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

    bool primary_finished() const {
      return primary_finished_;
    }

    void init() {
      if (cfg_->occupancy_images) {
        image_dir_ =
          std::filesystem::path(dpl::mpi::exec)
            .replace_filename("results")/cfg_->image.path.stem()/"images";
          // std::filesystem::path(dpl::mpi::exec).replace_filename("image")/settings_->image.path.stem();

        remove_all(image_dir_);
        create_directories(image_dir_);
      }

      state_.resize(pni_->connected_count(), pn_->throat_count());

      primary_.pc.reserve(500);
      primary_.kr.reserve(500);
      secondary_.pc.reserve(500);
      secondary_.kr.reserve(500);

      total_pore_volume_ = pni_->total_pore_volume(
        [this](voxel_t i) { return cfg_->image.darcy.info[img_->phase[i]].poro; });


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


    template <int cycle, int phase>
    double calc_relative(dpl::so_span<phase_t, const darcy_info> darcy, double theta) {
      dpl::so_uptr<phase_t, double> kr(darcy.size());

      for (auto* ptr = kr.data(); const auto& d : darcy)
        *ptr++ = d.kr[cycle][phase](d.pc_to_sw[cycle](state_.pc_global()));

      // for (phase_t i{0}; i < darcy.size(); ++i)
      //   kr[i] = darcy[i].kr[cycle][phase1](darcy[i].pc_to_sw[cycle](state_.pc_global()));

      filter_t<phase> filter{pni_, &img_->phase, &kr, &state_};
      term_multi_t term{
        std::bool_constant<phase>{},
        pni_,
        &state_,
        theta,
        [&](voxel_t i) {
          auto p = img_->phase[i];
          return kr[p]*darcy[p].perm;
        }
      };

      auto [nrows, mapping] = pni_->generate_mapping(cfg_->solver.decomposition, filter);
      auto [nvalues, input] = pni_->generate_pressure_input(nrows, std::move(mapping.forward), cfg_->macro_mult, term, filter);

      auto hash = pressure_cache::hash(nvalues, input);

      fmt::print("ph{} | hash: {:016x}", phase, hash);

      auto found = pressure_cache::cache().find(hash);

      if (!cfg_->solver.cache.use || found == pressure_cache::cache().end()) {
        dpl::so_uptr<net_t, HYPRE_Real> pressure(pni_->connected_count());
        auto decomposed_pressure = std::make_unique<HYPRE_Complex[]>(nrows);

        {
          using namespace std::chrono;

          // solve(input, {0, nrows - 1}, decomposed_pressure.get(), settings_->solver.tolerance, settings_->solver.max_iterations);

          auto t0 = high_resolution_clock::now();
          dpl::hypre::save_input(std::move(input), nrows, nvalues, mapping.block_rows, cfg_->solver.tolerance, cfg_->solver.max_iterations, 0, 0);
          auto t1 = high_resolution_clock::now();
          std::system(fmt::format("mpiexec -np {} \"{}\" -s", cfg_->solver.decomposition.prod(), dpl::mpi::exec).c_str()); // NOLINT(concurrency-mt-unsafe)
          auto t2 = high_resolution_clock::now();
          std::tie(decomposed_pressure, std::ignore, std::ignore) = dpl::hypre::load_values(nrows);

          std::cout << fmt::format(" | store: {} s, solve: {} s\n",
            duration_cast<seconds>(t1 - t0).count(), duration_cast<seconds>(t2 - t1).count());
        }
        
        for (HYPRE_BigInt i = 0; i < nrows; ++i)
          pressure[mapping.backward[i]] = decomposed_pressure[i];

        auto inlet = pni_->flow_rates(pressure, cfg_->macro_mult, term, filter).second;

        if (cfg_->solver.cache.save) {
          pressure_cache::cache()[hash] = inlet;
          pressure_cache::save();
        }

        return inlet;
      }

      std::cout << '\n';
      return found->second;
    }

    void launch_secondary(double absolute_rate, double theta) {
      using namespace std;
      using namespace attrib;
      using namespace dpl;
      using eq_tr = hydraulic_properties::equilateral_triangle;

      auto darcy_span = cfg_->image.darcy.span();
      auto darcy_r_cap = [](const darcy_info& d) { return 1/d.pc_to_sw[1].front().x; };

      so_uptr<net_t, bool> explored(pni_->connected_count());
      vector<bool> explored_throat(pn_->throat_count());

      displ_queue<false> queue;

      auto length_sum = 0.0;
      auto length_term = [this](auto i) { return volume(pn_, i)/eq_tr::area(r_ins(pn_, i)); };
      auto inv_vol0 = 0.0;

      auto has_films = [this](auto i) { return state_.config(i).layout() == phase_config::bulk_films(); };

      for (macro_t i{0}; i < pn_->node_count(); ++i)
        if (pni_->connected(i)) {
          if (has_films(pni_->net(i))) {
            length_sum += length_term(i);
            queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
          }
          else
            inv_vol0 += volume(pn_, i);
        }
      
      for (throat_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = adj(pn_, i); pni_->connected(l)) {
          if (has_films(i)) {
            length_sum += length_term(i);

            if (pn_->inlet() == r) {
              queue.insert(i, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, i)));
              explored_throat[i] = true;
            }
            else
              queue.insert(i, eq_tr::r_cap_snap_off(theta, r_ins(pn_, i)));
          }
          else
            inv_vol0 += volume(pn_, i);
        }

      so_uptr<phase_t, double> inv_poro(darcy_span.size(), 0);

      const double cell_volume = (pn_->physical_size/img_->dim()).prod();
      
      {
        idx3d_t ijk;
        auto& [i, j, k] = ijk;
        voxel_t idx1d{0};

        for (k = 0; k < img_->dim().z; ++k)
          for (j = 0; j < img_->dim().y; ++j)
            for (i = 0; i < img_->dim().x; ++i, ++idx1d)
              if (pni_->connected(idx1d)) {
                phase_t p = img_->phase[idx1d];
                
                if (has_films(pni_->net(idx1d))) {
                  inv_poro[p] += darcy_span[p].poro;
                  queue.insert(idx1d, darcy_r_cap(darcy_span[p]));  
                }
                else {
                  inv_vol0 += cell_volume*darcy_span[p].poro;      
                }
              }
      }


      
      auto eval_inv_volume = [&](double r_cap_global) {
        double vol = 0.0;
        for (phase_t i{0}; i < darcy_span.size(); ++i)
          vol += darcy_span[i].pc_to_sw[1](1/r_cap_global)*inv_poro[i];
        vol *= cell_volume;

        vol += inv_vol0 + length_sum*eq_tr::area_corners(theta, r_cap_global);

        return vol;
      };

      auto eval_pc_point = [&] {
        return vector2d{
          eval_inv_volume(state_.r_cap_global)/total_pore_volume_,
          state_.pc_global()
        };
      };

      using namespace dpl::graph;

      net_t outlet_idx = pni_->connected_count();
      et_ptr outlet_entry = get_entry(vertex_t(*outlet_idx), g_); 

      auto freeze_cluster = [&](et_ptr hdr) {
        auto area_corners = eq_tr::area_corners(theta, state_.r_cap_global);

        for (et_ptr et : range<et_traits>(hdr))
          if (is_loop_edge(et, g_)) {
            vertex_t v = get_vertex(et, g_);

            set_entry(v, nullptr, g_);

            net_t net{*v}; 

            if (has_films(net)) {  // NOLINT(clang-diagnostic-float-equal)
              state_.local(net) = state_.r_cap_global;

              if (pni_->is_macro(net)) {
                auto term = length_term(pni_->macro(net));
                length_sum -= term;
                inv_vol0 += area_corners*term;  
              }
              else {
                auto p = img_->phase[pni_->voxel(net)];
                auto poro = darcy_span[p].poro;

                inv_poro[p] -= poro;
                inv_vol0 += darcy_span[p].pc_to_sw[1](state_.pc_global())*cell_volume*poro;  
              }
            }

            if (pni_->is_macro(net)) {
              for (edge_t vw : g_.edges(v))
                if (pni_->is_macro(net_t{*target(vw, g_)})) {
                  if (throat_t t = de_to_throat_[vw]; has_films(t) && state_.mobile(t)) {
                    set_null_entry(vw, g_);
                    set_null_entry(opposite(vw, g_), g_);

                    state_.local(t) = state_.r_cap_global;
                    
                    auto term = length_term(t);
                    length_sum -= term;
                    inv_vol0 += area_corners*term;
                  }
                }
                else
                  set_null_entry(vw, g_);
            }
            else
              for (edge_t vu : g_.edges(v))
                set_null_entry(vu, g_);
          }
      };

      double last_kr_sw = 0.0;
      secondary_.kr.emplace_back(0, 0, 1);

      auto rel_calc_report = [&](double sw) {
        auto ph1_connected = [&] {
          for (auto [l, r] : pn_->throat_.span(adj))
            if (r == pn_->inlet() && pni_->connected(l)) // macro-inlet
              if (get_entry(vertex_t(*pni_->net(l)), g_))
                return true;

          {
            idx3d_t ijk;
            auto& [i, j, k] = ijk;
          
            for (k = 0; k < img_->dim().z; ++k)
              for (j = 0; j < img_->dim().y; ++j)
                if (voxel_t v{img_->idx_map(0, j, k)}; pni_->connected(v)) // darcy-inlet
                  if (get_entry(vertex_t(*pni_->net(v)), g_))
                    return true;
          }

          return false;
        };

        auto ph0_rate = calc_relative<1, 0>(darcy_span, theta);
        auto ph1_rate = ph1_connected() ? calc_relative<1, 1>(darcy_span, theta) : 0;
      
        secondary_.kr.emplace_back(sw, ph0_rate/absolute_rate, ph1_rate/absolute_rate);
      };

      auto write_occupancy_image = [&](double sw) {
        if (cfg_->occupancy_images)
          state_.write_occupancy_image<1>(
            image_dir_/fmt::format("secondary-{:.4f}.raw", sw), *pni_, theta, darcy_span);
      };

      if (darcy_span) {
        const double max_pc = primary_.pc.back().y;

        state_.r_cap_global = 1/max_pc;
        secondary_.pc.push_back(eval_pc_point());

        int steps = 10;

        const double topmost_pc = max(
          1/queue.front().r_cap,
          ranges::max(darcy_span | views::transform([=](const darcy_info& d) { return 1/darcy_r_cap(d); }))
        );

        if (const double step = pow(10, (log10(max_pc) - log10(topmost_pc))/steps); step > 0) {
          auto r_cap = state_.r_cap_global;
          
          for (int i = 0; i < steps; ++i) {
            r_cap *= step;
            if (auto p = eval_pc_point(); p.y < cfg_->max_pc) {
              state_.r_cap_global = r_cap;
              last_pc_point_ = p;
              secondary_.pc.push_back(last_pc_point_);
              if (i%2 != 0)
                write_occupancy_image(last_pc_point_.x);
            }
          }
        }

        ++progress_idx_;

        // if (true)
        { /* Relative permeabilty */

          set<double> pcs;
          for (const auto& d : darcy_span)
            for (auto pc : views::keys(d.pc_to_sw[1]))
              pcs.insert(pc);

          for (auto pc : pcs | std::views::reverse) {
            if (pc <= topmost_pc)
              break;

            if (auto sw = eval_inv_volume(1/pc)/total_pore_volume_;
              sw - last_kr_sw > cfg_->report.sw_of_kr &&
              pc < cfg_->max_pc)
            {
              state_.r_cap_global = 1/pc;
              last_kr_sw = sw;
              rel_calc_report(sw);

              ++progress_idx_;
            }
          }
        }
      }
      else {
        last_pc_point_ = eval_pc_point();
        secondary_.pc.push_back(last_pc_point_);
      }

      auto progress_percolation = [&]<typename Index>(Index i, double r_cap) {
        if (r_cap > state_.r_cap_global) {
          if (auto pc_point = eval_pc_point();
            abs(last_pc_point_.x - pc_point.x) > cfg_->report.sw_of_pc ||
            abs(log10(last_pc_point_.y) - log10(pc_point.y)) > pc_max_log_step_) {
            last_pc_point_ = pc_point;
            secondary_.pc.push_back(pc_point);
            write_occupancy_image(pc_point.x);
          }

          if (auto sw = eval_inv_volume(state_.r_cap_global)/total_pore_volume_;
            sw - last_kr_sw > cfg_->report.sw_of_kr)
          {
            last_kr_sw = sw;
            rel_calc_report(sw);
          }

          state_.r_cap_global = r_cap;
        }


        {
          state_.config(i) = phase_config{phase_config::phase0()};

          if constexpr (is_same_v<Index, net_t>) {
            vertex_t v{*i}; 

            for (edge_t vw : g_.edges(v))
              if (!is_null_entry(vw, g_) && !is_tree_edge(vw, g_))
                context_.non_tree_edge_remove(vw);

            for (edge_t vw : g_.edges(v))
              if (!is_null_entry(vw, g_) && !context_.tree_edge_split_and_reconnect(vw))  // remove remaining tree edges
                if (et_ptr hdr = et_algo::get_header(get_entry(target(vw, g_), g_));      // header of the target vertex
                  hdr != et_algo::get_header(outlet_entry))
                  freeze_cluster(hdr);

            set_entry(v, nullptr, g_);
          }
          else {
            if (edge_t vw = throat_to_de_[i]; !is_null_entry(vw, g_))
              if (is_tree_edge(vw, g_)) {
                if (!context_.tree_edge_split_and_reconnect(vw)) {
                  et_ptr hdr = et_algo::get_header(get_entry(target(vw, g_), g_));
                  if (hdr == et_algo::get_header(outlet_entry))
                    hdr = et_algo::get_header(get_entry(target(opposite(vw, g_), g_), g_));

                  freeze_cluster(hdr);
                }
              }
              else  // NOLINT(clang-diagnostic-dangling-else)
                context_.non_tree_edge_remove(vw);
          }
        }
      };


      // auto local_prog_idx = 0;

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
          net_t net = pni_->net(macro);

          if (state_.config(net).phase() == phase_config::phase0())
            continue;

          if (vertex_t v(*net); get_entry(v, g_)) {
            progress_percolation(net, r_cap);

            length_sum -= length_term(macro);
            inv_vol0 += volume(pn_, macro);

            for (edge_t vw : g_.edges(v)) {
              if (net_t w_net{*target(vw, g_)}; w_net == outlet_idx || pni_->is_macro(w_net)) { // macro-macro
                if (throat_t t = de_to_throat_[vw]; !explored_throat[t]) {
                  queue.insert(t, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, t)));
                  explored_throat[t] = true;
                }
              }
            }
          }
        }
        else if (elem == displ_elem::throat) {
          if (state_.config(local).phase() == phase_config::phase0())
            continue;

          auto [l, r] = adj(pn_, local);

          net_t l_net = pni_->net(l);

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

            length_sum -= length_term(local);
            inv_vol0 += volume(pn_, local);
          }
        }
        else { // voxel
          net_t net = pni_->net(voxel_t(local));  // NOLINT(CppTooWideScopeInitStatement)
          
          if (state_.r_cap_global == r_cap || get_entry(vertex_t(*net), g_))  // NOLINT(clang-diagnostic-float-equal)
          {
            progress_percolation(net, r_cap);
          }
        }
      }

      secondary_.pc.push_back(eval_pc_point());
      rel_calc_report(eval_inv_volume(state_.r_cap_global)/total_pore_volume_);

      ++progress_idx_;
    }

    void launch_primary(double absolute_rate, double theta) {
      // static constexpr bool to_sleep = false;

      pressure_cache::load();

      using namespace std;
      using namespace attrib;
      using namespace dpl;
      using eq_tr = hydraulic_properties::equilateral_triangle;

      auto darcy_span = cfg_->image.darcy.span();
      auto darcy_r_cap = [this, darcy_span](voxel_t i) { return 1/darcy_span[img_->phase[i]].pc_to_sw[0].front().x; };

      // if (darcy_span)
      pc_max_log_step_ = 0.25;
      // else {
      //   double min_r_cap_throat = numeric_limits<double>::max();
      //
      //   for (size_t i{0}; i < pn_->throat_count(); ++i)
      //     if (auto [l, r] = adj(pn_, i); pn_->inner_node(r) && pni_->connected(l))
      //       min_r_cap_throat = min(min_r_cap_throat, r_ins(pn_, i));
      //
      //   min_r_cap_throat = 0.95*eq_tr::r_cap_piston_with_films_valvatne(theta, min_r_cap_throat);
      //   pc_max_log_step_ = std::log10(0.075/min_r_cap_throat);
      // }

      so_uptr<net_t, bool> explored(pni_->connected_count());
      vector<bool> explored_throat(pn_->throat_count());

      displ_queue<true> queue;

      for (size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = adj(pn_, i); pn_->inlet() == r) // inlet macro-macro TODO: voxels inlet needed
          if (pni_->connected(l)) {
            queue.insert(i, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, i)));
            explored_throat[i] = true;
          }

      using namespace dpl::graph;

      net_t outlet_idx = pni_->connected_count();

      vector3d inv_vol_coefs{0};
      so_uptr<phase_t, double> inv_porosity(darcy_span.size(), 0);

      state_.r_cap_global = queue.front().r_cap;

      last_pc_point_ = {numeric_limits<double>::max(), numeric_limits<double>::min()};
      double last_kr_sw = 1.0;
      primary_.kr.emplace_back(1, 1, 0);

      const double cell_volume = (pn_->physical_size/img_->dim()).prod();
      auto eval_inv_volume = [&](double r_cap_global) {
        double vol = 0.0;

        for (phase_t i{0}; i < darcy_span.size(); ++i)
          vol += (1 - darcy_span[i].pc_to_sw[0](1/r_cap_global))*inv_porosity[i];
        vol *= cell_volume;

        vol += inv_vol_coefs[0] + inv_vol_coefs[2]*r_cap_global*r_cap_global;

        return vol;
      };

      auto eval_pc_point = [&] {
        return vector2d{
          1 - eval_inv_volume(state_.r_cap_global)/total_pore_volume_,
          state_.pc_global()
        };
      };

      auto rel_calc_report = [this, /*t0, */theta, absolute_rate, darcy_span](double sw) {
        // fmt::print("{}ms\n", duration_cast<milliseconds>(high_resolution_clock::now() - t0).count());
        primary_.kr.emplace_back(sw,
          calc_relative<0, 0>(darcy_span, theta)/absolute_rate,
          calc_relative<0, 1>(darcy_span, theta)/absolute_rate
          // sw,    
          // 1 - sw 
        ); // TODO: calculate when breakthrough
      };

      auto write_occupancy_image = [&](double sw) {
        if (cfg_->occupancy_images)
          state_.write_occupancy_image<0>(
            image_dir_/fmt::format("primary-{:.4f}.raw", sw), *pni_, theta, darcy_span);
      };

      auto progress_percolation = [&](auto idx, double r_cap) {
        if (r_cap < state_.r_cap_global) {
          auto pc_point = eval_pc_point();

          if (
            abs(last_pc_point_.x - pc_point.x) > cfg_->report.sw_of_pc ||
            abs(log10(last_pc_point_.y) - log10(pc_point.y)) > pc_max_log_step_)
          {
            last_pc_point_ = pc_point;
            primary_.pc.push_back(pc_point);
            write_occupancy_image(pc_point.x);
          }

          if (last_kr_sw - pc_point.x > cfg_->report.sw_of_kr) {
            last_kr_sw = pc_point.x;
            rel_calc_report(pc_point.x);
          }

          state_.r_cap_global = r_cap;
        }

        state_.config(idx) = phase_config::phase1_bulk_films();
      };

      


      progress_idx_ = 0;
      // getchar();

      so_uptr<phase_t, bool> slept{darcy_span.size()};
      for (/*progress_idx_ = 0*/; !queue.empty(); ++progress_idx_) {

        // {
        //   using namespace this_thread;
        //   using namespace chrono;
        //
        //   if (progress_idx_ < 500 && progress_idx_ % 25 == 0)
        //     sleep_for(milliseconds{250});
        //   // else if (progress_idx_ < 100)
        //   //   sleep_for(milliseconds{50});
        //   // else if (progress_idx_ < 300 /*1000*/)
        //   //   sleep_for(milliseconds{25});
        // }



        auto [elem, local, r_cap] = queue.front();

        if (1/r_cap > cfg_->max_pc) { /* terminating as max pc reached */
          last_pc_point_ = eval_pc_point();
          primary_.pc.push_back(last_pc_point_);
          write_occupancy_image(last_pc_point_.x);
          rel_calc_report(1 - eval_inv_volume(state_.r_cap_global)/total_pore_volume_);
          break;
        }

        queue.pop();

        if (elem == displ_elem::macro) {
          macro_t macro(local);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net = pni_->net(macro);

          // if (net == outlet_idx) {
          //   std::cout << "_____OUTLET______";
          // }

          progress_percolation(net, r_cap);

          inv_vol_coefs[0] += volume(pn_, macro);
          inv_vol_coefs[2] -= volume(pn_, macro)*eq_tr::area_corners(theta)/eq_tr::area(r_ins(pn_, macro));
          
          for (edge_t vu : g_.edges(vertex_t(*net)))
            if (net_t u_net_idx{*target(vu, g_)}; u_net_idx == outlet_idx || pni_->is_macro(u_net_idx)) { // macro-macro
              if (auto t_idx = de_to_throat_[vu]; !explored_throat[t_idx]) {
                queue.insert(t_idx, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, t_idx)));
                explored_throat[t_idx] = true;
              }
            }
            else { // macro-darcy
              if (!explored[u_net_idx]) {
                auto voxel = pni_->voxel(u_net_idx);
                queue.insert(voxel, darcy_r_cap(voxel));
                explored[u_net_idx] = true;
              }
            }
        }
        else if (elem == displ_elem::voxel) {
          voxel_t voxel(local);  // NOLINT(cppcoreguidelines-narrowing-conversions)
          auto net = pni_->net(voxel);

          darcy_invaded_ = true;

          {
            auto p = img_->phase[voxel];

            // if (to_sleep && !slept[p]) {
            //   slept[p] = true;
            //
            //   using namespace this_thread;
            //   using namespace chrono;
            //   
            //   sleep_for(milliseconds{1000});
            // }

            inv_porosity[p] += darcy_span[p].poro; // TODO -- merge with cell volume or remove poro from here
          }

          progress_percolation(net, r_cap); // TODO: phase_config::phase1_bulk_films(); of a voxel

          for (auto vu : g_.edges(vertex_t(*net))) {
            if (net_t u_net{*target(vu, g_)}; u_net != outlet_idx && !explored[u_net]) {
              if (pni_->is_macro(u_net)) { // macro
                auto u_macro = pni_->macro(u_net);
                queue.insert(u_macro, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, u_macro)));
                explored[u_net] = true;
              }
              else { // darcy
                auto u_voxel = pni_->voxel(u_net);
                queue.insert(u_voxel, darcy_r_cap(u_voxel));
                explored[u_net] = true;
              }
            }
          }
        }
        else { // throat
          auto [l, r] = adj(pn_, local);

          progress_percolation(local, r_cap);

          inv_vol_coefs[0] += volume(pn_, local);
          inv_vol_coefs[2] -= volume(pn_, local)*eq_tr::area_corners(theta)/eq_tr::area(r_ins(pn_, local));

          if (pn_->inner_node(r))
            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(r, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, r)));
              explored[r_net] = true;
            }

          if (auto l_net = pni_->net(l); !explored[l_net]) {
            queue.insert(l, eq_tr::r_cap_piston_with_films_valvatne(theta, r_ins(pn_, l)));
            explored[l_net] = true;
          }
        }
      }


      if (darcy_span) {
        auto last_pc = state_.pc_global();

        if (auto pc_point = eval_pc_point(); primary_.pc.back().x != pc_point.x) {  // NOLINT(clang-diagnostic-float-equal)
          primary_.pc.push_back(pc_point);
          write_occupancy_image(pc_point.x);
        }

        { /* Capillary pressure */
          constexpr auto steps = 10;

          using views::transform;

          auto darcy_max_pc = ranges::max(darcy_span | transform([](const darcy_info& d) { return d.pc_to_sw[0].back().x;  }));

          auto step = pow(10, (log10(darcy_max_pc) - log10(last_pc))/steps);

          for (auto i = 0; i < steps; ++i) {

            //  if (to_sleep &&/*i%2 == 0*/true) {
            //   ++progress_idx_;
            //   using namespace this_thread;
            //   using namespace chrono;
            //   sleep_for(milliseconds{1500});
            //   
            // }

            auto next_r_cap = state_.r_cap_global/step;

            if (1/next_r_cap > cfg_->max_pc) {
              state_.r_cap_global = 1/cfg_->max_pc;
              primary_.pc.push_back(eval_pc_point());
            
              if (i%2 != 0)
                write_occupancy_image(primary_.pc.back().x);

              break;
            }

            state_.r_cap_global = next_r_cap;
            primary_.pc.push_back(eval_pc_point());
          
            if (i%2 != 0) {
              write_occupancy_image(primary_.pc.back().x);
            }
          }

          ++progress_idx_;

          state_.r_cap_global = 1/last_pc;
        }

        
        if (true)
        { /* Relative permeabilty */

          set<double> pcs;
          for (const auto& d : darcy_span)
            for (auto pc : views::keys(d.pc_to_sw[0]))
              pcs.insert(pc);

          for (auto pc : pcs) {
            if (pc <= last_pc)
              continue;

            if (pc > cfg_->max_pc)
              break;

            if (auto total_sw = 1 - eval_inv_volume(1/pc)/total_pore_volume_;
              last_kr_sw - total_sw > cfg_->report.sw_of_kr)
            {
              state_.r_cap_global = 1/pc;
              last_kr_sw = total_sw;
              rel_calc_report(total_sw);

              ++progress_idx_;
            }
          }
        }
      }

      if (primary_.kr.back()[2] < 1)
        primary_.kr.emplace_back(0, 0, 1);

      ++progress_idx_;

      // getchar();

      primary_finished_ = true;

      {
        std::cout << '\n';
        generate_euler_tour();
        launch_secondary(absolute_rate, theta);
      }
    }
  };
}