#pragma once


#include "xpm/functions.h"

#include <dpl/units.hpp>
#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/qt/property_editor/PropertyItemsBase.hpp>
#include <dpl/qt/property_editor/QPropertyTreeView.hpp>
#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>
  
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






#undef LoadImage



namespace xpm
{
  


  class GlyphMapperFace
  {
  protected:
    vtkNew<vtkPoints> points_;
    std::vector<idx1d_t> original_indices_;

    vtkNew<vtkGlyph3DMapper> glyph_mapper_;
    vtkNew<vtkActor> actor_;
    vtkNew<vtkDoubleArray> color_arr_;

  public:
    auto* GetActor() const {
      return actor_.Get();
    }

    auto* GetGlyphMapper() const {
      return glyph_mapper_.Get();
    }

    auto& GetIndices() const {
      return original_indices_;
    }

    auto* GetColorArray() const {
      return color_arr_.Get();
    }
  };
  
  template<int face_idx>
  class GlyphMapperFaceGeneric : public GlyphMapperFace
  {
    using face = dpl::face_cubic<face_idx>;
    using dims = dpl::cdims<face::dim>;

    static vtkSmartPointer<vtkPolyData> Quad(double half_length = 0.5) {
      auto quad = vtkSmartPointer<vtkPolyData>::New();

      vtkNew<vtkPoints> points;
      quad->SetPoints(points);
      
      vtkNew<vtkCellArray> cells;
      quad->SetPolys(cells);
        
      v3d pos;
      pos[dims::e0] = 0;
      pos[dims::e1] = -half_length;
      pos[dims::e2] = -half_length;
      points->InsertNextPoint(pos);

      pos[dims::e1] = half_length;
      points->InsertNextPoint(pos);

      pos[dims::e2] = half_length;
      points->InsertNextPoint(pos);

      pos[dims::e1] = -half_length;
      points->InsertNextPoint(pos);

      if constexpr (face::is_upper) {
        vtkIdType indices[] = {0, 1, 2, 3};
        cells->InsertNextCell(4, indices);
      }
      else {
        vtkIdType indices[] = {3, 2, 1, 0};
        cells->InsertNextCell(4, indices);
      }

      return quad;
    }

    
  public:
    void Init(double half_length = 0.5) {
      
      // POLY_polydata__->SetPoints(POLY_points__);
      // POLY_polydata__->SetPolys(POLY_cells__);
      // POLY_polydata__->GetCellData()->SetScalars(color_arr_);




      // color_arr_->SetName("darcy_adj");

      vtkNew<vtkPolyData> polydata;
      polydata->GetPointData()->SetScalars(color_arr_);
      polydata->SetPoints(points_);

      glyph_mapper_->OrientOff();
      glyph_mapper_->SetScaleFactor(half_length/*1.0000*/); 
      glyph_mapper_->SetScaleModeToNoDataScaling();
      glyph_mapper_->SetInputData(polydata);
      glyph_mapper_->SetSourceData(Quad(/*half_length/2*/));
    }

    void Populate(const v3i& cells, const v3d& cell_size, const auto& filter) {
      idx3d_t map_idx{1, cells.x(), cells.x()*cells.y()};
      idx3d_t ijk;
      
      auto [e0, e1, e2] = dims::tie(ijk);
      auto [e0_count, e1_count, e2_count] = dims::tie(cells);
      
      auto adj_step = map_idx[dims::e0];

      v3d pos;
      
      idx1d_t idx1d;
      
      for (e2 = 0; e2 < e2_count; ++e2)
        for (e1 = 0; e1 < e1_count; ++e1) {
          e0 = 0;

          idx1d = ijk.dot(map_idx);

          if constexpr (!face::is_upper) {
            if (filter(idx1d)) {
              pos[dims::e0] = 0; //(0)*cell_size[e1_dim];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];

              points_->InsertNextPoint(pos);
              original_indices_.push_back(idx1d);
              // post(idx1d);

              // auto pts_count = POLY_points__->GetNumberOfPoints();
              //
              // pos[dims::e0] = 0;
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // pos[dims::e2] = (e2)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1 + 1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e2] = (e2 + 1)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);
              // vtkIdType indices[] = {pts_count + 3, pts_count + 2, pts_count + 1, pts_count + 0};
              // POLY_cells__->InsertNextCell(4, indices);

            }
          }
          
          for (; e0 < e0_count - 1; ++e0) {
            idx1d = ijk.dot(map_idx);
            bool filtered = filter(idx1d);

            auto adj_idx1d = idx1d + adj_step;
            bool adj_filtered = filter(adj_idx1d);

            if (filtered != adj_filtered) {
              pos[dims::e0] = (e0 + 1.0)*cell_size[dims::e0];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];



              // auto pts_count = POLY_points__->GetNumberOfPoints();
              //
              // pos[dims::e0] = (e0 + 1)*cell_size[dims::e0];
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // pos[dims::e2] = (e2)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1 + 1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e2] = (e2 + 1)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);

              if constexpr (face::is_upper) {
                if (filtered) {
                  points_->InsertNextPoint(pos);
                  original_indices_.push_back(idx1d);
                  // post(idx1d);

                  // vtkIdType indices[] = {pts_count + 0, pts_count + 1, pts_count + 2, pts_count + 3};
                  // POLY_cells__->InsertNextCell(4, indices);
                }
              }
              else {
                if (adj_filtered) {
                  points_->InsertNextPoint(pos);
                  original_indices_.push_back(adj_idx1d);
                  // post(adj_idx1d);


                  // vtkIdType indices[] = {pts_count + 3, pts_count + 2, pts_count + 1, pts_count + 0};
                  // POLY_cells__->InsertNextCell(4, indices);
                }
              }
            }
          }


          if constexpr (face::is_upper) {
            idx1d = ijk.dot(map_idx);
            
            if (filter(idx1d)) {
              pos[dims::e0] = (e0_count)*cell_size[dims::e0];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];
              
              points_->InsertNextPoint(pos);
              original_indices_.push_back(idx1d);
              // post(idx1d);





              // auto pts_count = POLY_points__->GetNumberOfPoints();
              //
              // pos[dims::e0] = (e0_count)*cell_size[dims::e0];
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // pos[dims::e2] = (e2)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1 + 1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e2] = (e2 + 1)*cell_size[dims::e2];
              // POLY_points__->InsertNextPoint(pos);
              //
              // pos[dims::e1] = (e1)*cell_size[dims::e1];
              // POLY_points__->InsertNextPoint(pos);
              // vtkIdType indices[] = {pts_count + 0, pts_count + 1, pts_count + 2, pts_count + 3};
              // POLY_cells__->InsertNextCell(4, indices);
            }
          }
        }

      color_arr_->SetNumberOfTuples(points_->GetNumberOfPoints());
    }
  };


  class ImageDataGlyphMapper
  {
  public:

    void Init(double half_length) {
      dpl::sfor<6>([=, this](auto i) {
        std::get<i>(faces_).Init(half_length);
      });
    }

    void Populate(const v3i& dims, const v3d& cell_size, const auto& filter) {
      auto start = std::chrono::high_resolution_clock::now();
      
      dpl::psfor<6>([=, this](auto i) {
        std::get<i>(faces_).Populate(dims, cell_size, filter);
        // std::cout << "\n\nFaces " << i;
      });

      auto stop = std::chrono::high_resolution_clock::now();
 
      cout << "\n\nFaces total time: " <<
        duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << endl;
    }

    std::tuple<
      GlyphMapperFaceGeneric<0>,
      GlyphMapperFaceGeneric<1>,
      GlyphMapperFaceGeneric<2>,
      GlyphMapperFaceGeneric<3>,
      GlyphMapperFaceGeneric<4>,
      GlyphMapperFaceGeneric<5>
    > faces_;
  };

  


  
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


    ImageDataGlyphMapper img_glyph_mapper_;








    image_data img_;

    


    bool use_cache = false;
    bool save_cache = true;

    

    


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
        auto color = QColor::fromHsl(coef*255, 175, 122);
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
          return std::get<0>(img_glyph_mapper_.faces_).GetActor()->GetProperty()->GetEdgeVisibility();
        };
        edges.set = [this](bool v) {
          dpl::sfor<6>([this, v](auto i) {
            static_cast<GlyphMapperFace&>(std::get<i>(img_glyph_mapper_.faces_)).GetActor()->GetProperty()->SetEdgeVisibility(v);
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




  public:
    // void LoadImage() {
    //   
    //
    //
    //   // auto image_path = 
    //   //   // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
    //   //   
    //   //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)"
    //   //   // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)"
    //   // ;
    //   //
    //   // auto velems_path = 
    //   //   R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
    //   //   // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
    //   // ;
    //
    //   
    //   // size_t image_size;
    //   
    // }


    
    
    void Init() {
      // auto image_path = R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Bmps252_6um.raw)";
      // auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)";
      // constexpr parse::image_dict input_spec{
      //   .solid = 1,       // dummy value, no '1' is in the image
      //   .pore = 255,
      //   .microporous = 0, // we read actual solid '0' as microporous
      // };
      

      auto image_path = R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)";
      auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)";
      constexpr parse::image_dict input_spec{
        .solid = 3,
        .pore = 0,
        .microporous = 2
      };


      //------------------------------

      v3i processors{1};

      if (auto proc_count = std::thread::hardware_concurrency();
        proc_count == 12)
        processors = {2, 2, 3};
      else if (proc_count == 32)
        processors = {4, 4, 2};



      
      static constexpr auto millidarcy = 1;
      static constexpr auto const_permeability = millidarcy*0.001*presets::darcy_to_m2;

      auto cache_path = std::format("cache/{}-pressure-{:.2f}mD.bin",
        std::filesystem::path(pnm_path).parent_path().filename().string(), const_permeability/presets::darcy_to_m2*1e3);


      InitGUI();

      InitLutNodeThroat(lut_node_throat_);
      InitLutPoreSolidMicro(lut_pore_solid_micro_);

      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      lut_pressure_->SetTableRange(0, 1);
      // lut_pressure_->SetNanColor(0.584314, 0.552624, 0.267419, 1);


      pore_network pn{pnm_path, pore_network::file_format::statoil};

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << "\n\nNetwork loaded";
        std::cout << (pn.eval_inlet_outlet_connectivity() ? " [CONNECTED]" : " [DISCONNECTED]");
      #endif

                

      pn.connectivity_flow_summary();



      {
        img_.read_image(image_path, input_spec);

        #ifdef XPM_DEBUG_OUTPUT
          std::cout << "\n\nImage phases read and array created";
        #endif

        img_.dim = std::round(std::cbrt(img_.size));
        img_.read_icl_velems(pnm_path);

        #ifdef XPM_DEBUG_OUTPUT
          std::cout << "\n\nVelems file read";
        #endif

        auto mapped_range = std::ranges::subrange{img_.velem.get(), img_.velem.get() + img_.size}
          | std::views::transform([](voxel_tag::velem x) { return *x; });

        InitLutVelems(lut_velem_, *std::ranges::max_element(mapped_range));

        #ifdef XPM_DEBUG_OUTPUT
          std::cout << "\n\nlut_velem_ created";
        #endif

        img_.eval_microporous_velem();

        #ifdef XPM_DEBUG_OUTPUT
          std::cout << "\n\nLoaded image";
        #endif
      }

     

      pore_network_image pni{pn, img_};

      #ifdef XPM_DEBUG_OUTPUT
      std::cout << "\n\nSTART_CONNECTIVITY ...";
      #endif

      pni.connectivity_inlet_outlet();

      #ifdef XPM_DEBUG_OUTPUT
      std::cout << "\n\nEND_CONNECTIVITY ///";
      #endif



      using namespace presets;

      std::vector<double> pressure;

      if (use_cache && std::filesystem::exists(cache_path)) {
        std::cout << "\n\nUsing CACHED pressure";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(double);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(nrows);
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(double));
      }
      else {
        std::cout << "\n\nparallel_partitioning_START";

        auto decomposition = pni.decompose_rows(processors);

        for (auto i = 0; i < processors.prod(); ++i)
          std::cout << std::format("\nblock {}, rows {}--{}, size {}",
            i,
            decomposition.rows_per_block[i].first,
            decomposition.rows_per_block[i].second,
            decomposition.rows_per_block[i].second - decomposition.rows_per_block[i].first + 1);

        std::cout << "\n\nparallel_partitioning_END";

        std::cout << "\n\nGeneratePressureInput START...";

        auto input = pni.generate_pressure_input(decomposition, const_permeability);

        std::cout << "\n\nGeneratePressureInput END|||";

        std::cout << "\n\nSave hypre input START...";

        using namespace boost::interprocess;
        
        {
          shared_memory_object smo{open_or_create, "xpm-hypre-input", read_write};
          dpl::hypre::save(input.get_ref(), input.nvalues, decomposition.rows_per_block, smo);
        }
        
        std::cout << "\n\nSave hypre input END|||";

        std::cout << "\n\nHypre solve MPI START...";
        
        auto start = std::chrono::high_resolution_clock::now();
        
        auto solve_result = std::system(
          std::format("mpiexec -n {} \"{}\" -s",
            processors.prod(), "xpm_project.exe"
          ).c_str());
        
        auto stop = std::chrono::high_resolution_clock::now();
        
        cout << "\n\nHypre solve MPI: " <<
          duration_cast<std::chrono::seconds>(stop - start).count() << "s END|||" << endl;
        
        shared_memory_object smo{open_only, "xpm-hypre-output", read_only};
        mapped_region mr{smo, read_only};
        
        auto opt_pressure = std::make_unique<double[]>(input.nrows);
        
        std::memcpy(opt_pressure.get(), mr.get_address(), input.nrows*sizeof(double));

        pressure.resize(input.nrows);
        for (auto i : dpl::range(input.nrows))
          pressure[decomposition.block_to_net[i]] = opt_pressure[i];

        std::cout << "\n\nPressure solved";


        std::filesystem::create_directory("cache");
        std::ofstream cache_stream(cache_path, std::ofstream::binary);
        cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(double)*pressure.size());
        std::cout << "\n\nPressure cached";
      }






      // std::vector<idx1d_t> _OLD_inner_parent(pn.node_count_);
      // std::iota(_OLD_inner_parent.begin(), _OLD_inner_parent.end(), 0);
      // std::vector<std::uint16_t> _OLD_inner_rank(pn.node_count_);
      // boost::disjoint_sets _OLD_inner_ds{_OLD_inner_rank.data(), _OLD_inner_parent.data()};
      // std::vector<bool> _OLD_inlet_connected(pn.node_count_);
      // std::vector<bool> _OLD_outlet_connected(pn.node_count_);
      //
      // auto _OLD_connected = [&](idx1d_t i) {
      //   auto rep_set = _OLD_inner_ds.find_set(i);
      //   return _OLD_inlet_connected[rep_set] && _OLD_outlet_connected[rep_set];
      // };
      //
      //
      //
      // {
      //   for (auto [l, r] : pn.throat_.range(attribs::adj))
      //     if (pn.inner_node(r))
      //       _OLD_inner_ds.union_set(l, r);
      //
      //
      //   for (auto [l, r] : pn.throat_.range(attribs::adj))
      //     if (r == pn.inlet())
      //       _OLD_inlet_connected[_OLD_inner_ds.find_set(l)] = true;
      //     else if (r == pn.outlet())
      //       _OLD_outlet_connected[_OLD_inner_ds.find_set(l)] = true;
      //
      //   idx1d_t disconnected_macro = 0;
      //   idx1d_t disconnected_macro_ALTER = 0;
      //   for (idx1d_t i = 0; i < macro_node_count; ++i) {
      //     if (!_OLD_connected(i))
      //       ++disconnected_macro;
      //
      //     if (!connected[i])
      //       ++disconnected_macro_ALTER;
      //   }
      //
      //   idx1d_t disconnected_darcy = 0;
      //   idx1d_t disconnected_darcy_ALTER = 0;
      //
      //   for (auto i : dpl::range(img_.size)) {
      //     if (img_.phase[i] == microporous && !_OLD_connected(macro_node_count + img_.darcy.index[i]))
      //       ++disconnected_darcy;
      //
      //     if (img_.phase[i] == microporous && !connected[macro_node_count + i])
      //       ++disconnected_darcy_ALTER;
      //   }
      //
      //   std::cout << std::format("\n\n Disconnected {} macro, {} darcy nodes", disconnected_macro, disconnected_darcy);
      //   std::cout << std::format("\n\n Disconnected {} macro, {} darcy nodes ALTER", disconnected_macro_ALTER, disconnected_darcy_ALTER);
      // }







      {
        // double inlet_flow_sum = 0;
        // double outlet_flow_sum = 0;
        //
        // for (idx1d_t k = 0; k < img_.dim.z(); ++k)
        //   for (idx1d_t j = 0; j < img_.dim.y(); ++j) {
        //     if (auto idx1d = map_idx(0, j, k);
        //       img_.phase[idx1d] == microporous)
        //       if (_OLD_connected(macro_node_count + img_.darcy.index[idx1d]))
        //         inlet_flow_sum += -2*cell_size.x()*const_permeability*(1 - pressure[macro_node_count + img_.darcy.index[idx1d]]);
        //
        //     if (auto idx1d = map_idx(img_.dim.x() - 1, j, k);
        //       img_.phase[idx1d] == microporous)
        //       if (_OLD_connected(macro_node_count + img_.darcy.index[idx1d]))
        //         outlet_flow_sum += -2*cell_size.x()*const_permeability*(pressure[macro_node_count + img_.darcy.index[idx1d]]);
        //   }
        //
        // for (idx1d_t i = 0; i < macro_throat_count; ++i) {
        //   auto [l, r] = pn.throat_[attribs::adj][i];
        //
        //   if (r == pn.inlet())
        //     if (_OLD_connected(l))
        //       inlet_flow_sum += coef_map(i)*(1 - pressure[l]);
        //   if (r == pn.outlet())
        //     if (_OLD_connected(l))
        //       outlet_flow_sum += coef_map(i)*(pressure[l]);
        // }
        //
        //
        //
        //
        //
        //
        // std::cout << std::format("\n\nMICROPOROUS_PERM={} mD\nINLET_PERM={} mD\nOUTLET_PERM={} mD\n",
        //   const_permeability/darcy_term*1000,
        //   -inlet_flow_sum/pn.physical_size.x()/darcy_term*1000,
        //   -outlet_flow_sum/pn.physical_size.x()/darcy_term*1000);
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        // for (auto [_, r] : pn.throat_.range(attribs::adj))
        //   if (!pn.inner_node(r))
        //     r -= pn.node_count_ - macro_node_count;
        //
        // pn.node_.resize(macro_node_count);
        // pn.node_count_ = macro_node_count;
        //
        // pn.throat_.resize(macro_node_count);
        // pn.throat_count_ = macro_throat_count;
      
      
      
      
      
        pni.flow_summary(pressure, const_permeability);
      }













      // std::ranges::fill(pressure, 0);









      
      

      

      

      

      


      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      assembly->AddPart(CreateNodeActor(pn, lut_pressure_, 
        [&](idx1d_t i) {
          return pni.connected_macro(i) ? pressure[pni.net_macro(i)] : std::numeric_limits<double>::quiet_NaN();
        }));

      assembly->AddPart(CreateThroatActor(pn, lut_pressure_, [&](size_t i) {
        auto [l, r] = pn.throat_[attribs::adj][i];

        return
          pni.connected_macro(l)
            ? pn.inner_node(r)
              ? pni.connected_macro(r) ? (pressure[pni.net_macro(l)] + pressure[pni.net_macro(r)])/2 : std::numeric_limits<double>::quiet_NaN()
              : r == pn.inlet() ? 1.0 : 0.0
            : std::numeric_limits<double>::quiet_NaN();
      }));

      renderer_->AddActor(assembly);
      
      std::cout << "\n\nNetwork actor created";


      


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

            img_glyph_mapper_.Populate(img_.dim, pn.physical_size/img_.dim, [&](idx1d_t idx1d) { return filter[idx1d]; });
              
            dpl::sfor<6>([&](auto face_idx) {
              GlyphMapperFace& face = std::get<face_idx>(img_glyph_mapper_.faces_);

              idx1d_t i = 0;
              for (auto idx1d : face.GetIndices()) {
                face.GetColorArray()->SetTypedComponent(i++, 0, 
                  // img_darcy_adj_arr[idx1d]
                  // velem_arr[idx1d].value
                  // total_ds.find_set(old_node_count + voxel_to_row_inc_map[idx1d]) == total_ds.find_set(total_parent.size() - 1) ? 0.5 : 1

                  img_.phase[idx1d] == microporous && pni.connected_darcy(idx1d)
                    ? pressure[pni.net_darcy(idx1d)]
                    : std::numeric_limits<double>::quiet_NaN()

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


      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());



      connect(
        qvtk_widget_,
        &QVTKWidgetRef::resized, this, [this]() {
          tidy_axes_.RefreshAxes();
        });
      
      
      
      // renderer_->GetActiveCamera()->AddObserver(vtkCommand::ModifiedEvent, this, &TidyAxes::RefreshAxes);
      
      tidy_axes_.SetFormat(".2e");

      double bounds[] = {
        0., pn.physical_size.x(),
        0., pn.physical_size.y(),
        0., pn.physical_size.z()};
      tidy_axes_.Build(bounds);
      // tidy_axes_.Build();
      // renderer_->ResetCamera();
    }
  };
}