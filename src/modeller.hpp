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

#include "invasion_task.hpp"

#include <list>

namespace xpm
{
  class modeller
  {
    pore_network pn_;
    image_data img_;
    pore_network_image pni_{pn_, img_};
    runtime_settings cfg_;


    struct
    {
      std::array<double, 2> perm_macro;   /* mD */
      std::array<double, 2> perm_total;   /* mD */

      double total_porosity;
      double effective_porosity;
    } petrophysics_summary_;

    double absolute_rate_;

    invasion_task invasion_task_{pni_, cfg_};

    template <bool primary>
    auto pc_to_string(fmt::format_string<const double&, const double&> f, auto sep, auto begin, auto end) const
    {
      using namespace std;

      list<string> rows;

      using dst = conditional_t<primary,
        front_insert_iterator<decltype(rows)>,
        back_insert_iterator<decltype(rows)>>;

      ranges::transform(
        primary ? invasion_task_.primary().pc : invasion_task_.secondary().pc,
        dst(rows), [f](const dpl::vector2d& p) { return fmt::format(f, p.x(), p.y()); });

      stringstream ss;
      ss << begin << rows.front();
      for (const auto& s : rows | views::drop(1))
        ss << sep << s;
      ss << end << '\n';

      return ss.str();
    }

    
    
  public:
    auto& petrophysics_summary() const {
      return petrophysics_summary_;
    }

    auto petrophysics_json() const {
      nlohmann::json j;

      j["poro"] = petrophysics_summary_.effective_porosity;

      {
        auto& perm = petrophysics_summary_.perm_total;
        j["perm"] = (perm[0] + perm[1])/2;
      }

      {
        auto& pc = j["cap_press"];

        {
          nlohmann::json primary;
          for (auto& p : invasion_task_.primary().pc | std::views::reverse)
            primary.push_back(p);
          pc.push_back(primary);
        }

        pc.push_back(invasion_task_.secondary().pc);
      }

      {
        auto& kr = j["rel_perm"];

        {
          nlohmann::json k0;
          nlohmann::json ph0;
          nlohmann::json ph1;

          for (const auto& p : invasion_task_.primary().kr | std::views::reverse) {
            ph0.push_back(dpl::vector2d{p.x(), p.y()});
            ph1.push_back(dpl::vector2d{p.x(), p.z()});
          }
                
          k0.push_back(ph0);
          k0.push_back(ph1);
          kr.push_back(k0);
        }

        {
          nlohmann::json k1;
          nlohmann::json ph0;
          nlohmann::json ph1;

          for (auto& p : invasion_task_.secondary().kr) {
            ph0.push_back(dpl::vector2d{p.x(), p.y()});
            ph1.push_back(dpl::vector2d{p.x(), p.z()});
          }

          k1.push_back(ph0);
          k1.push_back(ph1);
          kr.push_back(k1);
        }
      }

      return j;
    }

    template <bool primary>
    auto kr_to_string(fmt::format_string<const double&, const double&, const double&> f, auto sep, auto begin, auto end) const
    {
      using namespace std;

      list<string> rows;
        
      using dst = conditional_t<primary,
        front_insert_iterator<decltype(rows)>,
        back_insert_iterator<decltype(rows)>>;

      ranges::transform(
        primary ? invasion_task_.primary().kr : invasion_task_.secondary().kr,
        dst(rows), [f](const dpl::vector3d& p) { return fmt::format(f, p.x(), p.y(), p.z()); });

      stringstream ss;
      ss << begin << rows.front();
      for (const auto& s : rows | views::drop(1))
        ss << sep << s;
      ss << end;

      return ss.str();
    }

    template<bool primary>
    auto pc_to_plain() const {
      return modeller::pc_to_string<primary>("{:.6f}\t{:.6e}", "\n", "", "");
    }

    template<bool primary>
    auto pc_to_json() const {
      return modeller::pc_to_string<primary>("[{:.6f}, {:.6e}]", ", ", "[", "]");
    }

    template<bool primary>
    auto kr_to_plain() const {
      return modeller::kr_to_string<primary>("{:.6f}\t{:.6e}\t{:.6e}", "\n", "", "\n");
    }

    template<bool primary, bool phase1>
    auto kr_to_json() const {
      if constexpr (phase1)
        return modeller::kr_to_string<primary>("[{0:.6f}, {2:.6e}]", ", ", "[", "]");
      else
        return modeller::kr_to_string<primary>("[{0:.6f}, {1:.6e}]", ", ", "[", "]");
    }


    auto& pni() { return pni_; }
    auto& cfg() const { return cfg_; }
    auto absolute_rate() const { return absolute_rate_; }
    
    auto extract_network() {
      namespace fs = std::filesystem;

      auto filename = cfg_.image.pnextract_filename();

      if (auto copy_path = "pnextract"/filename; absolute(copy_path) != absolute(cfg_.image.path)) {
        // if (settings_.image.grey) {
        //   throw std::exception("not implemented (grey scale).");
        //
        //   // transform(
        //   //   settings_.image.path, settings_.image.size, 0,
        //   //   copy_path, settings_.image.size, [this](std::uint8_t v) {
        //   //     return
        //   //       v == 0 ? settings_.image.phases.pore :
        //   //       v == 255 ? settings_.image.phases.solid :
        //   //       settings_.image.phases.micro.get<std::uint8_t>();
        //   //   });
        //   //
        //   // settings_.darcy.read_grey(settings_.image.path, settings_.image.size);
        // }
        // else
          copy(cfg_.image.path, copy_path, fs::copy_options::update_existing);
      }

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      current_path(prev/"pnextract");
      auto network_dir = cfg_.image.path.stem();

      if (std::ranges::all_of(files, [&](std::string_view file) { return exists(network_dir/file); }))
        std::cout << "using cached network\n\n";
      else {
        std::cout << "=========== pnextract begin ===========\n" << std::flush;

        std::system( // NOLINT(concurrency-mt-unsafe)
          fmt::format("{} {}", fs::current_path()/"pnextract", filename).c_str());
      
        create_directory(network_dir);
        for (fs::path file : files)
          rename(file, network_dir/file);

        remove(fs::path{"_VElems.mhd"});

        std::cout << "=========== pnextract end ===========\n\n" << std::flush;
      }

      network_dir = absolute(network_dir);

      current_path(prev);

      return network_dir;
    }

    auto& get_invasion_task() {
      return invasion_task_;
    }

    void init(nlohmann::json& j) {
      cfg_.load(j);


      // for (voxel_ns::phase_t i{0}; i < cfg_.image.darcy.count; ++i) {
      //   for (int c = 0; c < 2; ++c)
      //     for (auto& p : cfg_.image.darcy.info[i].pc_to_sw[c] /*cfg_.image.darcy.info[i].pc_to_sw*/)
      //       for (int k = 0; k < *i; ++k)
      //         // p.x() *= 1.0;
      //         p.x() *= 1.2;    
      //
      //
      //
      //   // cfg_.image.darcy.info[i].kr = cfg_.primary.kr;
      // }


      std::filesystem::create_directories("cache");

      // std::locale::global(std::locale("en_US.UTF-8"));
    }

    void init(std::filesystem::path input) {
      if (!exists(input)) {
        std::cout << "Configuration file: ";
        std::cin >> input;

        if (!exists(input) && !input.has_extension())
          input.replace_extension("json");
          // std::cout << ".json";

        std::cout << '\n';
      }

      auto j = nlohmann::json::parse(std::ifstream{input}, nullptr, true, true);
      init(j.is_array() ? j[0] : j);
    }

    void prepare() {
      using namespace dpl;

      std::cout << fmt::format("image path: {}\n\n", cfg_.image.path);

      auto pnm_path = extract_network()/"";

      // auto begin_init_time = clock::now();

      //------------------------------

      // vector3i procs{1};
      //
      // if (settings_.solver.decomposition)
      //   procs = *settings_.solver.decomposition;
      // else {
      //   if (auto proc_count = std::thread::hardware_concurrency();
      //     proc_count == 12)
      //     procs = {2, 2, 3};
      //   else if (proc_count == 24)
      //     procs = {4, 3, 2};
      //   else if (proc_count == 32)
      //     procs = {4, 4, 2};
      //   else if (proc_count == 48)
      //     procs = {4, 4, 3};
      // }

      pn_.read_from_text_file(pnm_path);

      #ifdef XPM_DEBUG_OUTPUT
        std::cout
          << fmt::format("network\n  nodes: {:L}\n  throats: {:L}\n", pn_.node_count(), pn_.throat_count())
          << (pn_.eval_inlet_outlet_connectivity() ? "  (connected)" : "  (disconected)") << '\n';
      #endif

        
      petrophysics_summary_.perm_macro = 
      #ifdef _WIN32
        pn_.connectivity_flow_summary(cfg_.solver.tolerance, cfg_.solver.max_iterations, cfg_.macro_mult);
      #else
        pn_.connectivity_flow_summary_MPI(cfg_.solver.tolerance, cfg_.solver.max_iterations, cfg_.macro_mult);
      #endif

      std::cout << '\n';

      

      petrophysics_summary_.total_porosity = img_.read_image(cfg_.image);

      // petrophysics_summary_.porosity = img_.verify_darcy(settings_);

      img_.read_icl_velems(pnm_path);

      std::cout << "connectivity (isolated components)...";

      petrophysics_summary_.effective_porosity = pni_.evaluate_isolated(cfg_.image.darcy.info);
      
      std::cout << fmt::format(
        " done\n"
        "  eff. porosity: {:.2f}% \n"
        "  macro: {:L}\n"
        "  darcy: {:L}\n\n",
        petrophysics_summary_.effective_porosity*100,
        pni_.connected_macro_count(),
        *pni_.connected_count() - *pni_.connected_macro_count());
    }

    auto compute_pressure() {
      // using namespace dpl;
      //
      // std::cout << fmt::format("image path: {}\n\n", settings_.image.path);
      //
      // auto pnm_path = extract_network()/"";
      //
      // // auto begin_init_time = clock::now();
      //
      // //------------------------------
      //
      //
      // pn_.read_from_text_file(pnm_path);
      //
      // #ifdef XPM_DEBUG_OUTPUT
      //   std::cout
      //     << fmt::format("network\n  nodes: {:L}\n  throats: {:L}\n", pn_.node_count(), pn_.throat_count())
      //     << (pn_.eval_inlet_outlet_connectivity() ? "  (connected)" : "  (disconected)") << '\n';
      // #endif
      //
      // #ifdef _WIN32
      //   pn_.connectivity_flow_summary(settings_.solver.tolerance, settings_.solver.max_iterations);
      // #else
      //   pn_.connectivity_flow_summary_MPI(settings_.solver.tolerance, settings_.solver.max_iterations);
      // #endif
      //
      // std::cout << '\n';
      //
      // img_.dict = settings_.image.phases;
      //
      // img_.read_image("pnextract"/settings_.image.pnextract_filename());
      // img_.set_dim(settings_.loaded ? settings_.image.size : std::round(std::cbrt(*img_.size())));
      //
      // {
      //   using namespace voxel_prop;
      //
      //   strong_array<phase_t, bool> occurance;
      //   strong_array<phase_t, idx1d_t> count(0);
      //   for (auto p : img_.phase.span(img_.size())) {
      //     occurance[p] = true;
      //     ++count[p];
      //   }
      //
      //   std::list<phase_t> missing;
      //   for (int i = 0; i < 256; ++i)
      //     if (phase_t p(i);  // NOLINT(clang-diagnostic-implicit-int-conversion)
      //       settings_.image.phases.is_darcy(p) && occurance[p] && settings_.darcy.poro_perm[p].is_nan()) {
      //       missing.push_back(p);
      //     }
      //
      //   // { "value": 118, "porosity": 0.15, "permeability": 1.061493 }
      //
      //   if (!missing.empty()) {
      //     std::cout << "unassigned porosity and permeability: [\n";
      //
      //     auto print = [](auto iter) {
      //       fmt::print(R"(  {{ "value": {}, "porosity": null, "permeability": null }})", *iter);
      //     };
      //
      //     auto iter = missing.begin();
      //     auto end = missing.end();
      //
      //     print(iter);
      //     while (++iter != end) {
      //       std::cout << ",\n";
      //       print(iter);
      //     }
      //     std::cout << "\n]\n";
      //
      //     throw config_exception("missing microporosity values.");
      //   }
      // }
      //
      //
      // img_.read_icl_velems(pnm_path);
      //
      // std::cout << "connectivity (isolated components)...";
      //
      // pni_.evaluate_isolated();
      //
      // std::cout << fmt::format(" done\n  macro: {:L}\n  voxel: {:L}\n\n",
      //   pni_.connected_macro_count(),
      //   *pni_.connected_count() - *pni_.connected_macro_count());

      
      using namespace dpl;

      vector3i procs{1};
      
      if (cfg_.solver.decomposition)
        procs = *cfg_.solver.decomposition;
      else {
        if (auto proc_count = std::thread::hardware_concurrency();
          proc_count == 12)
          procs = {2, 2, 3};
        else if (proc_count == 24)
          procs = {4, 3, 2};
        else if (proc_count == 32)
          procs = {4, 4, 2};
        else if (proc_count == 48)
          procs = {4, 4, 3};
      }

      so_uptr<net_t, HYPRE_Real> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;




      std::cout << "matrix\n";
      fmt::print("  decompose {} {} procs...", procs.prod(), procs);

      auto [nrows, mapping] = pni_.generate_mapping(procs);

      std::cout
        << " done\n"
        << "  input build...";


      system::print_memory("PRE input", hypre::hypre_print_level);

      auto [nvalues, input] = pni_.generate_pressure_input(nrows, std::move(mapping.forward), cfg_.macro_mult, single_phase_conductance{&pn_,
        [this](voxel_t i) { return cfg_.image.darcy.info[img_.phase[i]].perm; }
      });

      system::print_memory("POST input", hypre::hypre_print_level);

      std::cout << " done\n";

      std::filesystem::create_directory("cache");

      // fmt::print("\n\nHASH: {:x}\n\n", pressure_cache::hash(nvalues, input));

      if (
        auto cache_path = fmt::format("cache/{}-pressure-{:x}.bin", cfg_.image.path.stem(), pressure_cache::hash(nvalues, input));
        cfg_.solver.cache.use && std::filesystem::exists(cache_path))
      {
        std::cout << "  using cache\n";
        pressure.resize(net_t(nrows));
        std::ifstream{cache_path, std::ifstream::binary}
          .read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));  // NOLINT(cppcoreguidelines-narrowing-conversions)
      }
      else {
        fmt::print("  input store [{:.0f} MB]...", units::megabyte{
          units::byte{
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex))
          }
        });


        save_input(std::move(input), nrows, nvalues, mapping.block_rows,
          cfg_.solver.tolerance, cfg_.solver.max_iterations, cfg_.solver.aggressive_levels, cfg_.solver.print_level);

        system::print_memory("PRE xpm MPI", hypre::hypre_print_level);


        std::cout
          << " done\n"
          // << "decomposition: (" << processors << fmt::format(") = {} procs\n\n", processors.prod())
          << "  solve hypre MPI...";

        using clock = std::chrono::high_resolution_clock;

        auto start = clock::now();

        std::system(  // NOLINT(concurrency-mt-unsafe)
          fmt::format("mpiexec -np {} \"{}\" -s", procs.prod(), mpi::exec).c_str()); 
        
        auto stop = clock::now();
        
        std::cout << " done " << duration_cast<std::chrono::seconds>(stop - start).count() << "s\n";

        {
          std::unique_ptr<HYPRE_Complex[]> decomposed_pressure;
          std::tie(decomposed_pressure, residual, iters) = hypre::load_values(nrows);

          pressure.resize(net_t(nrows));

          for (HYPRE_BigInt i = 0; i < nrows; ++i)
            pressure[mapping.backward[i]] = decomposed_pressure[i];
        }

        if (cfg_.solver.cache.save) {
          std::ofstream{cache_path, std::ofstream::binary}
            .write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)*nrows);
          std::cout << "  cached\n";
        }
      }

      // std::cout << '\n';

      {
        auto [inlet, outlet] = pni_.flow_rates(pressure, cfg_.macro_mult, single_phase_conductance{&pn_, 
          [this](voxel_t i) {
            return cfg_.image.darcy.info[img_.phase[i]].perm;
          }
        });

        using namespace presets;

        petrophysics_summary_.perm_total = {
          inlet*(pn_.physical_size.x()/(pn_.physical_size.y()*pn_.physical_size.z()))/darcy_to_m2*1000,
          outlet*(pn_.physical_size.x()/(pn_.physical_size.y()*pn_.physical_size.z()))/darcy_to_m2*1000
        };

        std::cout << fmt::format(
          // "microprs perm: {:.6f} mD\n"
          "  inlet perm: {:.6f} mD\n"
          "  outlet perm: {:.6f} mD\n"
          "  residual: {:.4g}, iterations: {}\n\n",
          // settings_.darcy.perm_single/darcy_to_m2*1000,
          petrophysics_summary_.perm_total[0],
          petrophysics_summary_.perm_total[1],
          residual, iters);

        absolute_rate_ = inlet;
      }

      return pressure;
    }
  };

  // struct saturation_map
  // {
  //   pore_network_image* pni;
  //   dpl::so_uptr<net_t, double>* data;
  //
  //   auto operator()(voxel_t voxel) const {
  //
  //     auto sw = 0.0;
  //
  //     // if (pni->connected(voxel))
  //     //   if (auto net = pni->net(voxel); state.config(net).phase() == phase_config::phase1())
  //     //     sw = 1.0 - solve(pc_inv, 1/state.r_cap(net), dpl::extrapolant::flat);
  //     //
  //     // face.GetColorArray()->SetTypedComponent(i++, 0, map_satur(sw));
  //
  //     // return pni->connected(i) ? (*data)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
  //   }
  //
  //   auto operator()(macro_t i) const {
  //     // return pni->connected(i) ? (*data)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
  //   }
  //
  //   auto operator()(throat_t i) const {
  //     auto [l, r] = attrib::adj(pni->pn(), i);
  //
  //     return pni->connected(l)
  //       ? ((*data)[pni->net(l)] + (*data)[pni->net(r)])/2
  //       : std::numeric_limits<double>::quiet_NaN();
  //   }
  //
  //   // auto operator()(voxel_t i) const {
  //   //   return pni->connected(i) ? (*pressure)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
  //   //
  //   //   pni().connected(voxel_t{idx1d})
  //   //         ?
  //   //           /*0*/
  //   //           pressure[pni().net(voxel_t{idx1d})]
  //   //           
  //   //           // (log10(model_.settings().darcy.perm(img().phase[voxel_t{idx1d}])) - log10(minPERM.perm))/(
  //   //           //   
  //   //           //   log10(maxPERM.perm) - log10(minPERM.perm))/1.25+0.1  /*/2+0.25*/
  //   //
  //   //         : std::numeric_limits<double>::quiet_NaN()
  //   // }
  // };
}
