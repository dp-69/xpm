#pragma once

#include "invasion_task.hpp"


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

      list<string> rows;

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

    auto kr_to_json(auto primary) const {
      return kr_to_string(primary, "[{0:.6f}, {2:.6e}]", ", ", "[", "]");
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
        if (settings_.image.grey) {
          transform(
            settings_.image.path, settings_.image.size, 0,
            copy_path, settings_.image.size, [this](std::uint8_t v) {
              return
                v == 0 ? settings_.image.phases.pore :
                v == 255 ? settings_.image.phases.solid :
                settings_.image.phases.micro.get<std::uint8_t>();
            });

          settings_.darcy.read_grey(settings_.image.path, settings_.image.size);
        }
        else
          copy(settings_.image.path, copy_path, fs::copy_options::update_existing);
      }

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      current_path(prev/"pnextract");
      auto network_dir = settings_.image.path.stem();

      if (std::ranges::all_of(files, [&](std::string_view file) { return exists(network_dir/file); }))
        std::cout << "using cached network\n\n";
      else {
        std::cout << "=========== pnextract's network extraction begin ===========\n";

        std::system( // NOLINT(concurrency-mt-unsafe)
          fmt::format("{} {}", fs::current_path()/"pnextract", filename).c_str());
      
        create_directory(network_dir);
        for (fs::path file : files)
          rename(file, network_dir/file);

        remove(fs::path{"_VElems.mhd"});

        std::cout << "=========== pnextract's network extraction end ===========\n\n";
      }

      network_dir = absolute(network_dir);

      current_path(prev);

      return network_dir;
    }

    auto& invasion_task() {
      return invasion_task_;
    }

    void init(const nlohmann::json& j) {
      settings_.load(j);

      std::filesystem::create_directories("cache");

      std::locale::global(std::locale("en_US.UTF-8"));
    }

    void init() {
      auto j = nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true);
      
      init(
        j.is_array() ? j[0] : j
      /*nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true)*/);

      // if (std::filesystem::exists("config.json"))
      //   settings_.load(nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true));
      //
      // std::filesystem::create_directories("cache");
      //
      // std::locale::global(std::locale("en_US.UTF-8"));
    }

    auto compute_pressure() {
      std::cout << fmt::format("image path: {}\n\n", settings_.image.path);

      auto pnm_path = extract_network()/"";

      // auto begin_init_time = clock::now();

      //------------------------------

      dpl::vector3i processors{1};

      if (settings_.solver.decomposition)
        processors = *settings_.solver.decomposition;
      else {
        if (auto proc_count = std::thread::hardware_concurrency();
          proc_count == 12)
          processors = {2, 2, 3};
        else if (proc_count == 24)
          processors = {4, 3, 2};
        else if (proc_count == 32)
          processors = {4, 4, 2};
        else if (proc_count == 48)
          processors = {4, 4, 3};
      }

      auto cache_path = fmt::format("cache/{}-pressure-{:.2f}mD.bin",
        settings_.image.path.stem(), settings_.darcy.perm_single/presets::darcy_to_m2*1e3);


      pn_.read_from_text_file(pnm_path);

      #ifdef XPM_DEBUG_OUTPUT
        std::cout
          << fmt::format("network\n  nodes: {:L}\n  throats: {:L}\n", pn_.node_count(), pn_.throat_count())
          << (pn_.eval_inlet_outlet_connectivity() ? "(connected)" : "(disconected)") << '\n';
      #endif

      #ifdef _WIN32
        pn_.connectivity_flow_summary(settings_.solver.tolerance, settings_.solver.max_iterations);
      #else
        pn_.connectivity_flow_summary_MPI(settings_.solver.tolerance, settings_.solver.max_iterations);
      #endif

      std::cout << '\n';

      {
        img_.read_image("pnextract"/settings_.image.pnextract_filename()/*settings_.image.path*/, settings_.image.phases);
        img_.set_dim(settings_.loaded ? settings_.image.size : std::round(std::cbrt(*img_.size())));

        img_.read_icl_velems(pnm_path);

        // auto mapped_range =
        //   std::span(img_.velem.data(), img_.size) 
        // | std::views::transform([](voxel_property::velem_t x) { return *x; });

        // InitLutVelems(lut_velem_, *std::ranges::max_element(mapped_range)); // TODO max is not valid, should check value validity

        img_.eval_microporous();
      }
      
      std::cout << "connectivity (isolated components)...";

      pni_.evaluate_isolated();

      std::cout << fmt::format(" done\n  macro: {:L}\n  voxel: {:L}\n\n",
        pni_.connected_macro_count(),
        pni_.connected_count() - pni_.connected_macro_count());

      

      dpl::strong_vector<net_t, double> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;


      if (settings_.solver.cache.use && std::filesystem::exists(cache_path)) {
        std::cout << "using cached pressure\n\n";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios::end);
        auto nrows = is.tellg()/sizeof(HYPRE_Complex);
        is.seekg(0, std::ios::beg);
        pressure.resize(net_t(nrows));
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));
      }
      else {
        std::cout << "decomposition...";

        auto [_, mapping] = pni_.generate_mapping(processors);

        std::cout
          << " done\n\n"
          << "input matrix build...";

        idx1d_t nrows = *pni_.connected_count();
        auto [nvalues, input] = pni_.generate_pressure_input(nrows, mapping.forward, single_phase_conductance{&pn_,
          [this](voxel_t i) { return settings_.darcy.perm(i); }
        });

        std::cout << " done\n\n";

        std::cout << 
          fmt::format("store input matrix [{} MB]...", (
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)))/1024/1024);

        dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_.solver.tolerance, settings_.solver.max_iterations);

        std::cout
          << " done\n\n"
          << "decomposition: (" << processors << fmt::format(") = {} procs\n\n", processors.prod())
          << "hypre MPI solve...";

        using clock = std::chrono::high_resolution_clock;

        auto start = clock::now();
        
        std::system(  // NOLINT(concurrency-mt-unsafe)
          fmt::format("mpiexec -np {} \"{}\" -s", processors.prod(), dpl::hypre::mpi::mpi_exec).c_str()); 
        
        auto stop = clock::now();
        
        std::cout << " done " << duration_cast<std::chrono::seconds>(stop - start).count() << "s\n\n";

        {
          std::unique_ptr<HYPRE_Complex[]> decomposed_pressure;
          std::tie(decomposed_pressure, residual, iters) = dpl::hypre::mpi::load_values(nrows);

          pressure.resize(net_t(nrows));

          for (HYPRE_BigInt i = 0; i < nrows; ++i)
            pressure[mapping.backward[i]] = decomposed_pressure[i];
        }

        std::cout << "pressure solved\n\n";

        if (settings_.solver.cache.save) {
          std::filesystem::create_directory("cache");
          std::ofstream cache_stream(cache_path, std::ofstream::binary);
          cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)**pni_.connected_count());
          std::cout << "pressure cached\n\n";
        }
      }

      {
        
        auto [inlet, outlet] = pni_.flow_rates(pressure, single_phase_conductance{&pn_, 
          [this](voxel_t i) { return settings_.darcy.perm(i); }
        });

        using namespace presets;

        std::cout << fmt::format(
          "microprs perm: {:.6f} mD\n"
          "  inlet perm: {:.6f} mD\n"
          "  outlet perm: {:.6f} mD\n"
          "  residual: {:.4g}, iterations: {}\n\n",
          settings_.darcy.perm_single/darcy_to_m2*1000,
          inlet/pn_.physical_size.x()/darcy_to_m2*1000,
          outlet/pn_.physical_size.x()/darcy_to_m2*1000,
          residual, iters);

        absolute_rate_ = inlet;
      }

      return pressure;
      

      // auto assembly = vtkSmartPointer<vtkAssembly>::New();
      //
      // {
      //   vtkSmartPointer<vtkActor> actor;
      //
      //   std::tie(actor, macro_colors) = CreateNodeActor(pn_, lut_pressure_, 
      //     [&](macro_t i) {
      //       return pni_.connected(i) ? /*0*/ pressure[pni_.net(i)] : std::numeric_limits<double>::quiet_NaN();
      //     });
      //
      //   assembly->AddPart(actor);
      //   
      //   std::tie(actor, throat_colors) = CreateThroatActor(pn_, lut_pressure_, [&](std::size_t i) {
      //     auto [l, r] = pn_.throat_[attrib::adj][i];
      //
      //     return pni_.connected(l)
      //       ? /*0*/ (pressure[pni_.net(l)] + pressure[pni_.net(r)])/2
      //       : std::numeric_limits<double>::quiet_NaN();
      //   });
      //
      //   assembly->AddPart(actor);
      // }
      //
      // renderer_->AddActor(assembly);
      //
      // // std::ranges::fill(pressure, 0);
      //
      // {
      //   auto scale_factor = /*1.0*/pn_.physical_size.x()/img_.dim().x(); // needed for vtk 8.2 floating point arithmetics
      //       
      //   img_glyph_mapper_.Init(scale_factor);
      //   {
      //     std::vector<bool> filter(*img_.size());
      //     for (voxel_t i{0}; i < img_.size(); ++i)
      //       filter[*i] = img_.phase[i] == microporous;
      //
      //     cout << "3D faces...";
      //
      //     auto t0 = clock::now();
      //
      //     img_glyph_mapper_.Populate(img_.dim(), pn_.physical_size/img_.dim(),
      //       [&](idx1d_t idx1d) { return filter[idx1d]; });
      //
      //     std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
      //         
      //     dpl::sfor<6>([&](auto face_idx) {
      //       dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
      //
      //       idx1d_t i = 0;
      //       for (auto idx1d : face.GetIndices())
      //         face.GetColorArray()->SetTypedComponent(i++, 0, 
      //           pni_.connected(voxel_t{idx1d})
      //             ? /*0*/ pressure[pni_.net(voxel_t{idx1d})]
      //             : std::numeric_limits<double>::quiet_NaN()
      //         );
      //
      //       auto* glyphs = face.GetGlyphMapper();
      //       auto* actor = face.GetActor();
      //
      //       // glyphs->SetLookupTable(lut_velem_);
      //       glyphs->SetLookupTable(lut_pressure_);
      //       // glyphs->SetLookupTable(lut_image_phase_);
      //         
      //       glyphs->SetColorModeToMapScalars();
      //       glyphs->UseLookupTableScalarRangeOn();
      //       // glyphs->SetScalarModeToUsePointData();
      //       // glyphs->SetScalarModeToUseCellData();
      //         
      //       actor->SetMapper(glyphs);
      //
      //       actor->GetProperty()->SetEdgeVisibility(false);
      //       actor->GetProperty()->SetEdgeColor(dpl::vector3d{0.25});
      //           
      //       actor->GetProperty()->SetAmbient(0.5);
      //       actor->GetProperty()->SetDiffuse(0.4);
      //       actor->GetProperty()->BackfaceCullingOn();
      //         
      //       renderer_->AddActor(actor);
      //     });
      //   }
      // }
      //
      //
      //
      // tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // // // tidy_axes_.SetFormat(".2e");
      //
      // bounds_ = {
      //   0., pn_.physical_size.x(),
      //   0., pn_.physical_size.y(),
      //   0., pn_.physical_size.z()};


    }


  };
}
