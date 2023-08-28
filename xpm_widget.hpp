#pragma once


#include "functions.h"

#include <dpl/units.hpp>
#include <dpl/hypre/mpi_module.hpp>
#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/qt/property_editor/PropertyItemsBase.hpp>
#include <dpl/qt/property_editor/QPropertyTreeView.hpp>
#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>
#include <dpl/vtk/ImageDataGlyphMapper.hpp>
  
// #include <QWidget>
#include <QMainWindow>
#include <QSplitter>


#if (VTK_MAJOR_VERSION == 8)
  #include <QVTKOpenGLWidget.h>
  #include <QSurfaceFormat>
#elif (VTK_MAJOR_VERSION == 9)
  #include <QVTKOpenGLNativeWidget.h>
#endif


#include <vtkAssembly.h>
#include <vtkCellData.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGlyph3DMapper.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>

#include <boost/pending/disjoint_sets.hpp>

#include <algorithm>
#include <future>

#include <nlohmann/json.hpp>



#undef LoadImage



namespace xpm
{
  


  

  


  
  #if (VTK_MAJOR_VERSION == 8)
    using QVTKWidgetRef = QVTKOpenGLWidget;
  #elif (VTK_MAJOR_VERSION == 9)
    using QVTKWidgetRef = QVTKOpenGLNativeWidget;
  #endif

  
  class XPMWidget : public QMainWindow
  {
  Q_OBJECT

    using QPropertyTreeView = dpl::qt::property_editor::QPropertyTreeView;
    
    QPropertyTreeView* tree_view_ = nullptr; 

    
    
    vtkNew<vtkRenderer> renderer_; 
    vtkNew<vtkGenericOpenGLRenderWindow> render_window_;
    
    // vtkRenderWindowInteractor* interactor_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;

    vtkNew<vtkLookupTable> lut_pressure_;
    vtkNew<vtkLookupTable> lut_continuous_;
    
    vtkNew<vtkLookupTable> lut_velem_;

    vtkNew<vtkLookupTable> lut_pore_solid_micro_;
    vtkNew<vtkLookupTable> lut_node_throat_;


    dpl::vtk::ImageDataGlyphMapper<idx1d_t> img_glyph_mapper_;








    image_data img_;
    std::array<double, 6> bounds_ = {0, 100, 0, 100, 0, 100};


    bool use_cache = true;
    bool save_cache = true;

    dpl::graph::dc_graph dc_graph_;
    dpl::graph::dc_context<dpl::graph::dc_properties> dc_context_;

    

    


    // static auto CreateNetworkAssembly(const pore_network_model& pnm, vtkLookupTable* lut) {
    //   auto net = vtkSmartPointer<vtkAssembly>::New();
    //   net->AddPart(CreateNodeActor(pnm, lut, [&](pnm_idx i) { return pnm.node_[attribs::r_ins][i]; }));
    //   net->AddPart(CreateThroatActor(pnm, lut, [&](pnm_idx i) { return pnm.throat_[attribs::r_ins][i]; }));
    //   return net;
    // }




    


    
    
    

    static void InitLutNodeThroat(vtkLookupTable* lut) {
      lut->IndexedLookupOn();
      lut->SetNumberOfTableValues(2);
      lut->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut->SetAnnotation(vtkVariant(0), "Node");
      lut->SetTableValue(1, 0.784314, 0.752624, 0.467419);
      lut->SetAnnotation(vtkVariant(1), "Throat");
    }

    static void InitLutPoreSolidMicro(vtkLookupTable* lut) {
      lut->IndexedLookupOn();  
      lut->SetNumberOfTableValues(3);
      lut->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut->SetAnnotation(vtkVariant(0), "Pore");
      lut->SetTableValue(1, 0.5, 0.5, 0.5);
      lut->SetAnnotation(vtkVariant(1), "Solid");
      lut->SetTableValue(2, 0.8, 0.8, 0.8);
      lut->SetAnnotation(vtkVariant(2), "Microporous");
    }

    static void InitLutVelems(vtkLookupTable* lut, std::int32_t max) {
      auto count = max - (-2)/**min*/ + 1;

      

      lut->SetNumberOfTableValues(count);
      
      for (int32_t i = 0; i < count; ++i) {
        auto coef = static_cast<double>(i*45%count)/count;
        auto color = QColor::fromHsl(coef*255, 175, 122);  // NOLINT(cppcoreguidelines-narrowing-conversions)
        lut->SetTableValue(i, color.redF(), color.greenF(), color.blueF());  // NOLINT(clang-diagnostic-double-promotion)
      }

      lut->SetTableValue(0, 0.7, 0.7, 0.7);     // -2 velem value
      lut->SetTableValue(1, 0.7, 0.7, 0.7);     // -1 velem value
      lut->SetTableRange(-2/**min*/, max);
    }


    

    void InitGUI() {
      qvtk_widget_ = new QVTKWidgetRef;

      #if (VTK_MAJOR_VERSION == 8)
        qvtk_widget_->SetRenderWindow(render_window_);
      #elif (VTK_MAJOR_VERSION == 9)
        qvtk_widget_->setRenderWindow(render_window_);
      #endif
      
      render_window_->AddRenderer(renderer_);
      // interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(v3d{1});


      {
        tree_view_ = new QPropertyTreeView;

        auto* model = tree_view_->model();
        // auto* cat_vis = model->AddCategory("Visualisation");

        dpl::qt::property_editor::ItemFunctor<bool> edges;

        edges.name = "Edges";
        edges.get = [this] {
          return static_cast<bool>(std::get<0>(img_glyph_mapper_.faces_).GetActor()->GetProperty()->GetEdgeVisibility());
        };
        edges.set = [this](bool v) {
          dpl::sfor<6>([this, v](auto i) {
            static_cast<dpl::vtk::GlyphMapperFace<idx1d_t>&>(std::get<i>(img_glyph_mapper_.faces_)).GetActor()->GetProperty()->SetEdgeVisibility(v);
            render_window_->Render();
          });
        };
        
        model->AddItem(/*cat_vis, */edges);
        
        auto* hsplit = new QSplitter{Qt::Horizontal};
        hsplit->addWidget(tree_view_);
        hsplit->addWidget(qvtk_widget_);
        hsplit->setStretchFactor(0, 0);
        hsplit->setStretchFactor(1, 1);
      
        setCentralWidget(hsplit);
      }

      // setCentralWidget(qvtk_widget_);
    }



    auto ProcessImage(const std::filesystem::path& p) {
      namespace fs = std::filesystem;

      auto filename = p.filename();

      if (auto copy_fn = "pnextract"/filename; fs::absolute(copy_fn) != fs::absolute(p))
        copy(p, copy_fn, fs::copy_options::update_existing);

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      fs::current_path(prev/"pnextract");
      auto network_dir = p.stem();

      if (std::ranges::all_of(files, [&](const auto& file) { return fs::exists(network_dir/file); })) {
        std::cout << "using cached network\n\n";
      }
      else {
        std::cout << "=========== pnextract's network extraction begin ===========\n";

        /*auto value = */std::system( // NOLINT(concurrency-mt-unsafe)
          fmt::format("{} {}", fs::current_path()/"pnextract", filename).c_str());
      
        fs::create_directory(network_dir);
        for (fs::path f : files)
          fs::rename(f, network_dir/f);

        fs::remove("_VElems.mhd");

        std::cout << "=========== pnextract's network extraction end ===========\n\n";
      }

      network_dir = fs::absolute(network_dir);

      fs::current_path(prev);

      return network_dir;
    }


  public:
    void Init() {
      std::locale::global(std::locale("en_US.UTF-8"));

      InitGUI();


      ComputePressure();


      renderer_->ResetCamera(bounds_.data());
      renderer_->GetActiveCamera()->Zoom(0.70);
      // renderer_->ResetCameraClippingRange();
      
      tidy_axes_.Init(renderer_.Get());

      connect(qvtk_widget_, &QVTKWidgetRef::resized, this, [this]() { tidy_axes_.RefreshAxes(); });

      // tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // tidy_axes_.SetFormat(".2e");
      tidy_axes_.Build(bounds_.data());
    }

    void ComputePressure() {
      using clock = std::chrono::high_resolution_clock;
      using seconds = std::chrono::seconds;

      startup_settings startup;

      auto has_config = std::filesystem::exists("config.json");
      if (has_config)
        startup.load(nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true));


      std::cout << fmt::format("image path: {}\n\n", startup.image.path);

      auto pnm_path = ProcessImage(startup.image.path)/"";

      auto begin_init_time = clock::now();

      //------------------------------

      v3i processors{1};

      if (startup.solver.decomposition)
        processors = *startup.solver.decomposition;
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

      auto const_permeability = startup.microporous_perm*0.001*presets::darcy_to_m2;

      
      auto cache_path = fmt::format("cache/{}-pressure-{:.2f}mD.bin",
        startup.image.path.stem(), const_permeability/presets::darcy_to_m2*1e3);



      InitLutNodeThroat(lut_node_throat_);
      InitLutPoreSolidMicro(lut_pore_solid_micro_);

      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      lut_pressure_->SetTableRange(0, 1);
      // lut_pressure_->SetNanColor(0.584314, 0.552624, 0.267419, 1);


      pore_network pn{pnm_path, pore_network::file_format::statoil};

      #ifdef XPM_DEBUG_OUTPUT
        std::cout
          << "network\n  "
          << (pn.eval_inlet_outlet_connectivity() ? "(connected)" : "(disconected)") << '\n';
      #endif


      #ifdef _WIN32
        pn.connectivity_flow_summary(startup.solver.tolerance, startup.solver.max_iterations);
      #else
        pn.connectivity_flow_summary_MPI(startup.solver.tolerance, startup.solver.max_iterations);
      #endif

      std::cout << '\n';

      {
        img_.read_image(startup.image.path, startup.image.phases);
        img_.dim = has_config ? startup.image.size : std::round(std::cbrt(img_.size));

        img_.read_icl_velems(pnm_path);

        auto mapped_range =
          std::ranges::subrange{img_.velem.get(), img_.velem.get() + img_.size}
        | std::views::transform([](voxel_tag::velem x) { return *x; });

        InitLutVelems(lut_velem_, *std::ranges::max_element(mapped_range));

        img_.eval_microporous();

        // #ifdef XPM_DEBUG_OUTPUT
        //   std::cout << "image loaded\n\n";
        // #endif
      }

      pore_network_image pni{pn, img_};

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << "connectivity (isolated components)...";
      #endif

      pni.evaluate_isolated();

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << " done\n\n";
      #endif



      using namespace presets;

      std::vector<HYPRE_Complex> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;


      if (use_cache && std::filesystem::exists(cache_path)) {
        std::cout << "using cached pressure\n\n";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(HYPRE_Complex);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(nrows);
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));
      }
      else {
        std::cout << "decomposition...";

        auto decomposition = pni.decompose_rows(processors);

        // for (auto i = 0; i < processors.prod(); ++i)
        //   std::cout << std::format("\nblock {}, rows {}--{}, size {}",
        //     i,
        //     decomposition.rows_per_block[i].first,
        //     decomposition.rows_per_block[i].second,
        //     decomposition.rows_per_block[i].second - decomposition.rows_per_block[i].first + 1);

        std::cout
          << " done\n\n"
          << "input matrix build...";

        auto [nrows, nvalues, input] = pni.generate_pressure_input(decomposition, const_permeability);

        std::cout
          << " done\n\n"
          // << "pre hypre time: " <<
          // duration_cast<std::chrono::milliseconds>(clock::now() - begin_init_time).count() << "ms\n\n"
        ;

        std::cout << 
          fmt::format("save input matrix [{} MB]...", (
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)))/1024/1024);

        dpl::hypre::mpi::save(input, nrows, nvalues, decomposition.blocks, startup.solver.tolerance, startup.solver.max_iterations);

        std::cout
          << " done\n\n"
          << "decomposition: (" << processors << fmt::format(") = {} procs\n\n", processors.prod())
          << "hypre MPI solve...";
        
        auto start = clock::now();
        
        std::system(  // NOLINT(concurrency-mt-unsafe)
          fmt::format("mpiexec -np {} \"{}\" -s", processors.prod(), dpl::hypre::mpi::mpi_exec).c_str()); 
        
        auto stop = clock::now();
        
        cout << " done " << duration_cast<std::chrono::seconds>(stop - start).count() << "s\n\n";





        // auto [decomposed_pressure, out_residual, out_iters] = dpl::hypre::mpi::load_values(nrows);
        // residual = out_residual;
        // iters = out_iters;

        {
          std::unique_ptr<HYPRE_Complex[]> decomposed_pressure;
          std::tie(decomposed_pressure, residual, iters) = dpl::hypre::mpi::load_values(nrows);

          pressure.resize(nrows);
          for (auto i : dpl::range(nrows))
            pressure[decomposition.decomposed_to_net[i]] = decomposed_pressure[i];
        }

        std::cout << "pressure solved\n\n";

        if (save_cache) {
          std::filesystem::create_directory("cache");
          std::ofstream cache_stream(cache_path, std::ofstream::binary);
          cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)*pressure.size());
          std::cout << "pressure cached\n\n";
        }
      }


      pni.flow_summary(pressure.data(), const_permeability);
      std::cout << fmt::format("  residual: {:.4g}, iterations: {}\n\n", residual, iters);

      {
        using namespace dpl::graph;
        
        auto dc_props = dc_properties{dc_graph_};

        auto t0 = clock::now();
        std::cout << "graph...";
        dc_graph_ = pni.generate_dc_graph();
        std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
          duration_cast<seconds>(clock::now() - t0).count(),
          dc_graph_.vertex_count(),
          dc_graph_.edge_count());

        auto t1 = clock::now();
        std::cout << "euler tour..."; // dfs for et-only, fast etnte
        dc_context_.init_with_dfs(dc_graph_, dc_props);
        std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t1).count());
      }






      // std::ranges::fill(pressure, 0);

      

      


      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      // assembly->AddPart(CreateNodeActor(pn, lut_pressure_, 
      //   [&](macro_idx i) {
      //     return
      //       et_algo::get_header(dc_props.get_entry(dc_graph_.get_vertex(pn.node_count() + 1))) ==
      //       et_algo::get_header(dc_props.get_entry(dc_graph_.get_vertex(*i)))
      //         ? 0
      //         : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
      //
      //     // return pni.connected(i) ? pressure[pni.net(i)] : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
      //   }));
      //
      // assembly->AddPart(CreateThroatActor(pn, lut_pressure_, [&](size_t i) {
      //   auto [l, r] = pn.throat_[attribs::adj][i];
      //
      //   return
      //     pni.connected(l)
      //       ? pn.inner_node(r)
      //         ? pni.connected(r) ? (pressure[pni.net(l)] + pressure[pni.net(r)])/2 : std::numeric_limits<HYPRE_Complex>::quiet_NaN()
      //         : r == pn.inlet() ? 1.0 : 0.0
      //       : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
      // }));

      assembly->AddPart(CreateNodeActor(pn, lut_pressure_, 
        [&](macro_idx i) {
          return pni.connected(i) ? pressure[pni.net(i)] : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
        }));
      
      assembly->AddPart(CreateThroatActor(pn, lut_pressure_, [&](size_t i) {
        auto [l, r] = pn.throat_[attribs::adj][i];
      
        return
          pni.connected(l)
            ? pn.inner_node(r)
              ? pni.connected(r) ? (pressure[pni.net(l)] + pressure[pni.net(r)])/2 : std::numeric_limits<HYPRE_Complex>::quiet_NaN()
              : r == pn.inlet() ? 1.0 : 0.0
            : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
      }));

      renderer_->AddActor(assembly);
      
      // std::cout << "\n\nNetwork actor created";

      {
        
          
        {
          auto scale_factor = /*1.0*/pn.physical_size.x()/img_.dim.x(); // needed for vtk 8.2 floating point arithmetics
            
          img_glyph_mapper_.Init(scale_factor);
          {
            std::vector<bool> filter(img_.size);
            for (idx1d_t i = 0; i < img_.size; ++i) {
              // auto rep_set = inner_ds.find_set(macro_node_count + img_.darcy.index[i]);
              filter[i] = img_.phase[i] == microporous
              // && !con_io[pn.node_count_ + i]
              // && (con_io[pn.node_count_ + i] != connected(macro_node_count + img_.darcy.index[i]))
              // && inlet_connected[rep_set] && outlet_connected[rep_set]
              // && total_ds.find_set(old_node_count + voxel_to_row_inc_map[i]) == outlet_set
              ;
            }

            cout << "3D faces...";

            auto t0 = clock::now();

            img_glyph_mapper_.Populate(img_.dim, pn.physical_size/img_.dim, [&](idx1d_t idx1d) { return filter[idx1d]; });

            std::cout << fmt::format(" done {}s\n", duration_cast<seconds>(clock::now() - t0).count());
              
            dpl::sfor<6>([&](auto face_idx) {
              dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);

              idx1d_t i = 0;
              for (auto idx1d : face.GetIndices()) {
                voxel_idx v_idx{idx1d};

                face.GetColorArray()->SetTypedComponent(i++, 0, 
                  // img_darcy_adj_arr[idx1d]
                  // velem_arr[idx1d].value
                  // total_ds.find_set(old_node_count + voxel_to_row_inc_map[idx1d]) == total_ds.find_set(total_parent.size() - 1) ? 0.5 : 1

                  img_.phase[idx1d] == microporous && pni.connected(v_idx)
                    ? pressure[pni.net(v_idx)]
                    : std::numeric_limits<HYPRE_Complex>::quiet_NaN()

                  // phase_arr[idx1d].value
                );
              }

              auto* glyphs = face.GetGlyphMapper();
              auto* actor = face.GetActor();

              // glyphs->SetLookupTable(lut_velem_);
              glyphs->SetLookupTable(lut_pressure_);
              // glyphs->SetLookupTable(lut_image_phase_);
              
              glyphs->SetColorModeToMapScalars();
              glyphs->UseLookupTableScalarRangeOn();
              // glyphs->SetScalarModeToUsePointData();
              // glyphs->SetScalarModeToUseCellData();
              
              actor->SetMapper(glyphs);

              actor->GetProperty()->SetEdgeVisibility(/*false*/false);
              actor->GetProperty()->SetEdgeColor(v3d{0.25} /*0, 0, 0*/);
                
              actor->GetProperty()->SetAmbient(0.5);
              actor->GetProperty()->SetDiffuse(0.4);
              actor->GetProperty()->BackfaceCullingOn();
              
              renderer_->AddActor(actor);
            });
          }
        }
      }


      // renderer_->ResetCamera();
      //
      // tidy_axes_.Init(renderer_.Get());
      //
      // connect(
      //   qvtk_widget_,
      //   &QVTKWidgetRef::resized, this, [this]() {
      //     tidy_axes_.RefreshAxes();
      //   });
      //
      // // renderer_->GetActiveCamera()->AddObserver(vtkCommand::ModifiedEvent, this, &TidyAxes::RefreshAxes);
      //
      tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // // tidy_axes_.SetFormat(".2e");

      bounds_ = {
        0., pn.physical_size.x(),
        0., pn.physical_size.y(),
        0., pn.physical_size.z()};

      // tidy_axes_.Build(bounds);
      // // tidy_axes_.Build();
      // // renderer_->ResetCamera();
    }
  };
}



// std::filesystem::path image_path = R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Bmps-v0s255_252x252x252_6p0um.raw)";
      // parse::image_dict input_spec{
      //   .pore = 0,
      //   .solid = 1,       // dummy value, no '1' is in the image
      //   .microporous = 255, // we read actual solid as microporous
      // };

      // std::filesystem::path image_path = R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Est-v0m2s3_500x500x500_4p0um.raw)";
      // constexpr parse::image_dict input_spec{
      //   .solid = 3,
      //   .pore = 0,
      //   .microporous = 2
      // };
     
      // HYPRE_Real tolerance = 1.e-9; HYPRE_Int max_iterations = 1000;