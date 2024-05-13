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
    runtime_settings settings_;
    double absolute_rate_;
    invasion_task invasion_task_{pni_, settings_};


    auto pc_to_string(auto primary, fmt::format_string<const double&, const double&> f, auto sep, auto begin, auto end) const {
      using namespace std;

      std::list<string> rows;

      using dst = conditional_t<primary,
        front_insert_iterator<decltype(rows)>,
        back_insert_iterator<decltype(rows)>>;

      ranges::transform(
        primary ? invasion_task_.primary().pc : invasion_task_.secondary().pc,
        dst(rows), [f](const dpl::vector2d& p) { return format(f, p.x(), p.y()); });

      stringstream ss;
      ss << begin << rows.front();
      for (const auto& s : rows | views::drop(1))
        ss << sep << s;
      ss << end << '\n';

      return ss.str();
    }

    
    
  public:
    auto kr_to_string(
      auto primary,
      fmt::format_string<const double&, const double&, const double&> f, auto sep, auto begin, auto end) const {
      using namespace std;
      list<string> rows;
        
      using dst = conditional_t<primary,
        front_insert_iterator<decltype(rows)>,
        back_insert_iterator<decltype(rows)>>;

      ranges::transform(
        primary ? invasion_task_.primary().kr : invasion_task_.secondary().kr,
        dst(rows), [f](const dpl::vector3d& p) { return format(f, p.x(), p.y(), p.z()); });

      stringstream ss;
      ss << begin << rows.front();
      for (const auto& s : rows | views::drop(1))
        ss << sep << s;
      ss << end;

      return ss.str();
    }

    auto pc_to_plain(auto primary) const {
      return pc_to_string(primary, "{:.6f}\t{:.6e}", "\n", "", "");
    }

    auto pc_to_json(auto primary) const {
      return pc_to_string(primary, "[{:.6f}, {:.6e}]", ", ", "[", "]");
    }

    auto kr_to_plain(auto primary) const {
      return kr_to_string(primary, "{:.6f}\t{:.6e}\t{:.6e}", "\n", "", "\n");
    }

    auto kr_to_json(auto primary, auto phase1) const {
      if constexpr (phase1)
        return kr_to_string(primary, "[{0:.6f}, {2:.6e}]", ", ", "[", "]");
      else
        return kr_to_string(primary, "[{0:.6f}, {1:.6e}]", ", ", "[", "]");
    }


    // auto& pn() { return pn_; }
    // auto& img() { return img_; }
    auto& pni() { return pni_; }
    const auto& settings() { return settings_; }
    auto absolute_rate() { return absolute_rate_; }
    
    auto extract_network() {
      namespace fs = std::filesystem;

      auto filename = settings_.image.pnextract_filename();

      if (auto copy_path = "pnextract"/filename; absolute(copy_path) != absolute(settings_.image.path)) {
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
          copy(settings_.image.path, copy_path, fs::copy_options::update_existing);
      }

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      current_path(prev/"pnextract");
      auto network_dir = settings_.image.path.stem();

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

    void init(const nlohmann::json& j) {
      settings_.load(j);

      std::filesystem::create_directories("cache");

      std::locale::global(std::locale("en_US.UTF-8"));
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

      std::cout << fmt::format("image path: {}\n\n", settings_.image.path);

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

      #ifdef _WIN32
        pn_.connectivity_flow_summary(settings_.solver.tolerance, settings_.solver.max_iterations);
      #else
        pn_.connectivity_flow_summary_MPI(settings_.solver.tolerance, settings_.solver.max_iterations);
      #endif

      std::cout << '\n';

      img_.dict = settings_.image.phases;

      img_.read_image("pnextract"/settings_.image.pnextract_filename());
      img_.set_dim(settings_.loaded ? settings_.image.size : std::round(std::cbrt(*img_.size())));

      {
        using namespace voxel_prop;

        strong_array<phase_t, bool> occurances;
        strong_array<phase_t, idx1d_t> count(0);
        for (auto p : img_.phase.span(img_.size())) {
          occurances[p] = true;
          ++count[p];
        }

        std::list<phase_t> missing;
        for (int i = 0; i < 256; ++i)
          if (phase_t p(i);  // NOLINT(clang-diagnostic-implicit-int-conversion)
            settings_.image.phases.is_darcy(p) && occurances[p] && settings_.darcy.poro_perm[p].is_nan()) {
            missing.push_back(p);
          }

        // { "value": 118, "porosity": 0.15, "permeability": 1.061493 }

        if (!missing.empty()) {
          std::cout << "unassigned porosity and permeability: [\n";

          auto print = [](auto iter) {
            fmt::print(R"(  {{ "value": {}, "porosity": null, "permeability": null }})", *iter);
          };

          auto iter = missing.begin();
          auto end = missing.end();

          print(iter);
          while (++iter != end) {
            std::cout << ",\n";
            print(iter);
          }
          std::cout << "\n]\n";

          throw config_exception("missing microporosity values.");
        }
      }


      img_.read_icl_velems(pnm_path);

      std::cout << "connectivity (isolated components)...";

      pni_.evaluate_isolated();

      std::cout << fmt::format(" done\n  macro: {:L}\n  voxel: {:L}\n\n",
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
      
      if (settings_.solver.decomposition)
        procs = *settings_.solver.decomposition;
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

      strong_vector<net_t, HYPRE_Real> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;




      std::cout << "matrix\n";
      fmt::print("  decompose {} {} procs...", procs.prod(), procs);

      auto [nrows, mapping] = pni_.generate_mapping(procs);

      std::cout
        << " done\n"
        << "  input build...";


      system::print_memory("PRE input", hypre::hypre_print_level);

      auto [nvalues, input] = pni_.generate_pressure_input(nrows, std::move(mapping.forward), single_phase_conductance{&pn_,
        [this](voxel_t i) { return settings_.darcy.perm(img_.phase[i]); }
      });

      system::print_memory("POST input", hypre::hypre_print_level);

      std::cout << " done\n";

      std::filesystem::create_directory("cache");

      if (
        auto cache_path = fmt::format("cache/{}-pressure-{:x}.bin", settings_.image.path.stem(), pressure_cache::hash(nvalues, input));
        settings_.solver.cache.use && std::filesystem::exists(cache_path))
      {
        std::cout << "  using cache\n";
        pressure.resize(net_t(nrows));
        std::ifstream{cache_path, std::ifstream::binary}
          .read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));  // NOLINT(cppcoreguidelines-narrowing-conversions)
      }
      else {
        fmt::print("  input store [{} MB]...", units::megabyte{
          units::byte{
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex))
          }
        });


        save_input(std::move(input), nrows, nvalues, mapping.block_rows,
          settings_.solver.tolerance, settings_.solver.max_iterations, settings_.solver.aggressive_levels);

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

        if (settings_.solver.cache.save) {
          std::ofstream{cache_path, std::ofstream::binary}
            .write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)*nrows);
          std::cout << "  cached\n";
        }
      }

      // std::cout << '\n';

      {
        auto [inlet, outlet] = pni_.flow_rates(pressure, single_phase_conductance{&pn_, 
          [this](voxel_t i) { return settings_.darcy.perm(img_.phase[i]); }
        });

        using namespace presets;

        std::cout << fmt::format(
          // "microprs perm: {:.6f} mD\n"
          "  inlet perm: {:.6f} mD\n"
          "  outlet perm: {:.6f} mD\n"
          "  residual: {:.4g}, iterations: {}\n\n",
          // settings_.darcy.perm_single/darcy_to_m2*1000,
          inlet*(pn_.physical_size.x()/(pn_.physical_size.y()*pn_.physical_size.z()))/darcy_to_m2*1000,
          outlet*(pn_.physical_size.x()/(pn_.physical_size.y()*pn_.physical_size.z()))/darcy_to_m2*1000,
          residual, iters);

        absolute_rate_ = inlet;
      }

      return pressure;
    }
  };

  struct saturation_map
  {
    pore_network_image* pni;
    dpl::strong_vector<net_t, double>* data;

    auto operator()(voxel_t voxel) const {

      auto sw = 0.0;

      // if (pni->connected(voxel))
      //   if (auto net = pni->net(voxel); state.config(net).phase() == phase_config::phase1())
      //     sw = 1.0 - solve(pc_inv, 1/state.r_cap(net), dpl::extrapolant::flat);
      //
      // face.GetColorArray()->SetTypedComponent(i++, 0, map_satur(sw));

      // return pni->connected(i) ? (*data)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
    }

    auto operator()(macro_t i) const {
      // return pni->connected(i) ? (*data)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
    }

    auto operator()(throat_t i) const {
      auto [l, r] = attrib::adj(pni->pn(), i);

      return pni->connected(l)
        ? ((*data)[pni->net(l)] + (*data)[pni->net(r)])/2
        : std::numeric_limits<double>::quiet_NaN();
    }

    // auto operator()(voxel_t i) const {
    //   return pni->connected(i) ? (*pressure)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
    //
    //   pni().connected(voxel_t{idx1d})
    //         ?
    //           /*0*/
    //           pressure[pni().net(voxel_t{idx1d})]
    //           
    //           // (log10(model_.settings().darcy.perm(img().phase[voxel_t{idx1d}])) - log10(minPERM.perm))/(
    //           //   
    //           //   log10(maxPERM.perm) - log10(minPERM.perm))/1.25+0.1  /*/2+0.25*/
    //
    //         : std::numeric_limits<double>::quiet_NaN()
    // }
  };
}
