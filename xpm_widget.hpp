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

#include <algorithm>
#include <future>




#undef LoadImage



namespace xpm
{
  namespace parse
  {
    struct image_values
    {
      std::uint8_t solid;
      std::uint8_t pore;
      std::uint8_t microporous;
    };

    std::pair<std::unique_ptr<voxel_tag::phase[]>, std::size_t> read_image(const auto& image_path, image_values input_config) {
      std::ifstream is(image_path);
      is.seekg(0, std::ios_base::end);
      size_t image_size = is.tellg();
      is.seekg(0, std::ios_base::beg);
      auto phase_arr = std::make_unique<voxel_tag::phase[]>(image_size);
      is.read(reinterpret_cast<char*>(phase_arr.get()), image_size);

      size_t pore_voxels = 0;
      size_t solid_voxels = 0;
      size_t microporous_voxels = 0;

      for (size_t i = 0; i < image_size; ++i)
        if (auto& val = phase_arr[i];
          val.value == input_config.pore) {
          val = presets::pore;
          ++pore_voxels;
        }
        else if (val.value == input_config.solid) {
          val = presets::solid;
          ++solid_voxels;
        }
        else if (val.value == input_config.microporous) {
          val = presets::microporous;
          ++microporous_voxels;
        }

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << std::format("\n\ntotal: {}; pore: {}; solid: {}; microporous: {} voxels", image_size, pore_voxels, solid_voxels, microporous_voxels);
      #endif

      return {std::move(phase_arr), image_size};
    }
  }


  class GlyphMapperFace
  {
  protected:
    vtkNew<vtkPoints> points_;
    std::vector<pnm_idx> original_indices_;

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
      pnm_3idx map_idx{1, cells.x(), cells.x()*cells.y()};
      pnm_3idx ijk;
      
      auto [e0, e1, e2] = dims::tie(ijk);
      auto [e0_count, e1_count, e2_count] = dims::tie(cells);
      
      auto adj_step = map_idx[dims::e0];

      v3d pos;
      
      pnm_idx idx1d;
      
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
    
    std::unique_ptr<QPropertyTreeView> tree_view_factory_;
    QPropertyTreeView* tree_view_ = nullptr; 

    
    
    vtkNew<vtkRenderer> renderer_; 
    vtkNew<vtkGenericOpenGLRenderWindow> render_window_;
    
    vtkRenderWindowInteractor* interactor_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;

    vtkNew<vtkLookupTable> lut_pressure_;
    vtkNew<vtkLookupTable> lut_continuous_;
    
    vtkNew<vtkLookupTable> lut_velem_;

    vtkNew<vtkLookupTable> lut_pore_solid_micro_;
    vtkNew<vtkLookupTable> lut_node_throat_;



    // vtkNew<vtkActor> image_actor_;




    ImageDataGlyphMapper img_mapper;




    std::unique_ptr<voxel_tag::phase[]> phase_arr;
    std::unique_ptr<voxel_tag::velem[]> velem_arr;

    /*
     * relevant for microporous phase only. 
     *
     * connected pore node of a microporous voxel. Gives a cluster number [0, n].
     * -2 for pore voxels,
     * -1 for solid voxels not connected
     * 
     */
    std::unique_ptr<std::int32_t[]> adj_pore_of_darcy_arr;
    


    bool use_cache = true;
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


    static bool InletOutletConnectivity(pore_network_model& pnm) {
      std::vector<pnm_idx> parent(pnm.node_count_ + 2);
      std::vector<std::uint16_t> rank(pnm.node_count_ + 2);

      std::iota(parent.begin(), parent.end(), 0);

      boost::disjoint_sets ds{rank.data(), parent.data()};


      for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
        ds.union_set(n0, n1);

      return ds.find_set(pnm.inlet()) == ds.find_set(pnm.outlet());
    }

    void InitGUI() {
      qvtk_widget_ = new QVTKWidgetRef;

      #if (VTK_MAJOR_VERSION == 8)
        qvtk_widget_->SetRenderWindow(render_window_);
      #elif (VTK_MAJOR_VERSION == 9)
        qvtk_widget_->setRenderWindow(render_window_);
      #endif
      
      render_window_->AddRenderer(renderer_);
      interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(v3d{1});


      {
        tree_view_ = new QPropertyTreeView;

        auto* model = tree_view_->model();
        // auto* cat_vis = model->AddCategory("Visualisation");

        dpl::qt::property_editor::ItemFunctor<bool> edges;

        edges.name = "Edges";
        edges.get = [this] {
          return std::get<0>(img_mapper.faces_).GetActor()->GetProperty()->GetEdgeVisibility();
        };
        edges.set = [this](bool v) {
          dpl::sfor<6>([this, v](auto i) {
            static_cast<GlyphMapperFace&>(std::get<i>(img_mapper.faces_)).GetActor()->GetProperty()->SetEdgeVisibility(v);
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
    void LoadImage() {
      // auto image_path = R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)";
      // auto velems_path = R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)";
      // constexpr struct
      // {
      //   std::uint8_t solid = 1;       // dummy value, no '1' is in the image
      //   std::uint8_t pore = 255;
      //   std::uint8_t microporous = 0; // we read actual solid '0' as microporous
      // } input_spec;


      auto image_path = R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)";
      auto velems_path = R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)";
      constexpr parse::image_values input_spec{
        .solid = 3,
        .pore = 0,
        .microporous = 2
      };


      







              // auto image_path = 
              //   // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
              //   
              //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)"
              //   // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)"
              // ;
              //
              // auto velems_path = 
              //   R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
              //   // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
              // ;




      




      

      


      
      
   

      



      size_t image_size;
      std::tie(phase_arr, image_size) = parse::read_image(image_path, input_spec);

      std::cout << "\n\nImage phases read and array created";

      pnm_3idx dim = std::round(std::cbrt(image_size));

      {
        velem_arr = pore_network_model::read_icl_velems(velems_path, dim);

        std::cout << "\n\nVelems file read";

        auto mapped_range = std::ranges::subrange{velem_arr.get(), velem_arr.get() + image_size}
          | std::views::transform([](voxel_tag::velem x) { return x.value; });
        auto max = *std::ranges::max_element(mapped_range);
        auto count = max - (-2)/**min*/ + 1;

        std::cout << "\n\nRANGES MINMAX VELEMS";
        
        lut_velem_->SetNumberOfTableValues(count);

        for (int32_t i = 0; i < count; ++i) {
          auto coef = static_cast<double>(i*45%count)/count;
          auto color = QColor::fromHsl(coef*255, 175, 122);
          lut_velem_->SetTableValue(i, color.redF(), color.greenF(), color.blueF());  // NOLINT(clang-diagnostic-double-promotion)
        }

        lut_velem_->SetTableValue(0, 0.7, 0.7, 0.7);     // -2 velem value
        lut_velem_->SetTableValue(1, 0.7, 0.7, 0.7);     // -1 velem value
        lut_velem_->SetTableRange(-2/**min*/, max);
        
        std::cout << "\n\nlut_velem_ created";
        
        {
          adj_pore_of_darcy_arr = std::make_unique<std::int32_t[]>(dim.prod());

          pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};

          pnm_idx idx1d = 0;

          pnm_3idx ijk;
          auto& [i, j, k] = ijk;

          for (k = 0; k < dim.z(); ++k)
            for (j = 0; j < dim.y(); ++j)
              for (i = 0; i < dim.x(); ++i, ++idx1d) {
                int32_t adj = -1;
                
                if (phase_arr[idx1d] == presets::microporous)
                  dpl::sfor<3>([&](auto d) {
                    if (ijk[d] > 0)
                      if (auto adj_velem = velem_arr[idx1d - map_idx[d]].value; adj_velem >= 0)
                        adj = adj_velem;
          
                    if (ijk[d] < dim[d] - 1)
                      if (auto adj_velem = velem_arr[idx1d + map_idx[d]].value; adj_velem >= 0)
                        adj = adj_velem;
                  });

                adj_pore_of_darcy_arr[idx1d] = adj;
              }
        
          std::cout << "\n\nVelems_adj array produced and filled";
        }
      }
    }


    
    
    void Init() {
      // auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)";
      // v3i dim = 252;

      auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)";
      v3i dim = 500;

      v3i processors{1};

      if (auto proc_count = std::thread::hardware_concurrency();
        proc_count == 12)
        processors = {2, 2, 3};
      else if (proc_count == 32)
        processors = {4, 4, 2};


      using std::filesystem::path;


      static constexpr auto darcy_term = 9.869233e-13;
      static constexpr auto millidarcy = 1;
      static constexpr auto const_permeability = millidarcy*0.001*darcy_term;

      auto cache_path = std::format("cache/{}-pressure-{:.2f}mD.bin",
        path(pnm_path).parent_path().filename().string(), const_permeability/darcy_term*1e3);


      InitGUI();

      InitLutNodeThroat(lut_node_throat_);
      InitLutPoreSolidMicro(lut_pore_solid_micro_);

      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      lut_pressure_->SetTableRange(0, 1);


      pore_network_model pnm{pnm_path, pore_network_model::file_format::statoil};

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << "\n\nNetwork loaded";
        std::cout << (InletOutletConnectivity(pnm) ? " [CONNECTED]" : " [DISCONNECTED]");
      #endif

      LoadImage();

      #ifdef XPM_DEBUG_OUTPUT
        std::cout << "\n\nLoaded image";
      #endif


      auto cell_size = pnm.physical_size/dim;
      pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};

      auto old_node_count = pnm.node_count_;
      auto old_throat_count = pnm.throat_count_;

      auto coef_map = [&](size_t i) {
        using eq_tri = geometric_properties::equilateral_triangle_properties;
        using namespace attribs;

        auto [n0, n1] = pnm.throat_[adj][i];

        if (!pnm.inner_node(n0) && pnm.inner_node(n1) && n1 >= old_node_count) {
          auto coef = -2*cell_size.x()*const_permeability;
          return coef;
        }

        if (!pnm.inner_node(n1) && pnm.inner_node(n0) && n0 >= old_node_count) {
          auto coef = -2*cell_size.x()*const_permeability;
          return coef;
        }

        if (pnm.inner_node(n0) && pnm.inner_node(n1) && n0 >= old_node_count && n1 >= old_node_count) {
          auto coef = -cell_size.x()*const_permeability;
          return coef;
        }

        auto coef = -1.0/(
          (pnm.inner_node(n0) ? pnm.throat_[length0][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n0])) : 0.0) +
          (pnm.inner_node(n1) ? pnm.throat_[length1][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n1])) : 0.0) +
          pnm.throat_[length][i]/eq_tri::conductance(eq_tri::area(pnm.throat_[r_ins][i])));

        return coef;
      };




      std::vector<double> pressure;
      auto voxel_to_row_inc_map = std::make_unique<pnm_idx[]>(dim.prod());


      {
        auto image_size = dim.prod<pnm_idx>();

        pnm_idx darcy_outer_count = 0; // darcy nodes connected to macroscopic pore nodes
        pnm_idx darcy_inner_count = 0; // darcy nodes not connected to marcoscopic pore nodes

        for (pnm_idx i = 0; i < image_size; ++i)
          if (phase_arr[i] == presets::microporous) {
            if (adj_pore_of_darcy_arr[i] >= 0)
              ++darcy_outer_count;
            else
              ++darcy_inner_count;
          }


        pnm_idx darcy_darcy_throats = 0;



        // pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};
        //
        // pnm_idx idx1d = 0;
        //
        // pnm_3idx ijk;
        // auto& [i, j, k] = ijk;
        //
        // for (k = 0; k < dim.z(); ++k)
        //   for (j = 0; j < dim.y(); ++j)
        //     for (i = 0; i < dim.x(); ++i, ++idx1d) {
        //       int32_t adj = -1;
        //       
        //       if (phase_arr[idx1d] == presets::microporous)
        //         dpl::sfor<3>([&](auto d) {
        //           if (ijk[d] > 0)
        //             if (auto adj_velem = velem_arr[idx1d - map_idx[d]].value; adj_velem >= 0)
        //               adj = adj_velem;
        //
        //           if (ijk[d] < dim[d] - 1)
        //             if (auto adj_velem = velem_arr[idx1d + map_idx[d]].value; adj_velem >= 0)
        //               adj = adj_velem;
        //         });
        //
        //       adj_pore_of_darcy_arr[idx1d] = adj;
        //     }


        {
          pnm_idx idx1d = 0;
          pnm_3idx ijk;
          auto& [i, j, k] = ijk;

          for (k = 0; k < dim.z(); ++k)
            for (j = 0; j < dim.y(); ++j)
              for (i = 0; i < dim.x(); ++i, ++idx1d)
                if (phase_arr[idx1d] == presets::microporous)
                  dpl::sfor<3>([&](auto d) {
                    if (ijk[d] < dim[d] - 1)
                      if (phase_arr[idx1d + map_idx[d]] == presets::microporous)
                        ++darcy_darcy_throats;
                  });
        }

        auto darcy_inlet_outlet = 0;

        for (pnm_idx k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j) {
            if (phase_arr[map_idx.dot({0, j, k})] == presets::microporous)
              ++darcy_inlet_outlet;
            if (phase_arr[map_idx.dot({dim.x() - 1, j, k})] == presets::microporous)
              ++darcy_inlet_outlet;
          }

        std::cout << "\n\nPRE_RESIZE";
          
        pnm.node_.resize(pnm.node_count_ + darcy_outer_count + darcy_inner_count);
        pnm.throat_.resize(pnm.throat_count_ + darcy_outer_count + darcy_darcy_throats + darcy_inlet_outlet);

        std::cout << "\n\nPOST_RESIZE";
          
        pnm_idx new_throat_incr = 0;
          
        for (auto& [n0, n1] : pnm.throat_.range(attribs::adj)) {
          if (n0 == pnm.inlet() || n0 == pnm.outlet())
            n0 += darcy_outer_count + darcy_inner_count;

          if (n1 == pnm.inlet() || n1 == pnm.outlet())
            n1 += darcy_outer_count + darcy_inner_count;
        }

          
        pnm_idx voxel_to_row_inc = 0;
          
        for (pnm_idx idx1d = 0, k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j)
            for (pnm_idx i = 0; i < dim.x(); ++i, ++idx1d) {
              if (phase_arr[idx1d] == presets::microporous) {
                auto adj_pore_idx = adj_pore_of_darcy_arr[idx1d];

                if (adj_pore_idx >= 0) {
                  auto new_node_idx = old_node_count + (voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++);
                    
                  pnm.node_[attribs::r_ins][new_node_idx] = cell_size.x()/2;
                  pnm.node_[attribs::pos][new_node_idx] = cell_size*(v3d{i, j, k} + 0.5);
                    
                  auto new_throat_idx = pnm.throat_count_ + new_throat_incr++;
                    
                  pnm.throat_[attribs::adj][new_throat_idx] = {new_node_idx, adj_pore_idx};
                  pnm.throat_[attribs::r_ins][new_throat_idx] = cell_size.x()/4;//min_r_ins;
                    
                  pnm.throat_[attribs::length0][new_throat_idx] = 0;
                  pnm.throat_[attribs::length][new_throat_idx] =
                    (pnm.node_[attribs::pos][new_node_idx] - pnm.node_[attribs::pos][adj_pore_idx]).length();
                  pnm.throat_[attribs::length1][new_throat_idx] = 0;
                }
                else {
                  auto new_node_idx = old_node_count + (voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++);
                  
                  pnm.node_[attribs::r_ins][new_node_idx] = cell_size.x()/2;
                  pnm.node_[attribs::pos][new_node_idx] = cell_size*(v3d{i, j, k} + 0.5);
                }
              }
            }


        std::cout << "\n\nNODES_VALUES";


        {
          pnm_idx idx1d = 0;
          pnm_3idx ijk;
          auto& [i, j, k] = ijk;

          for (k = 0; k < dim.z(); ++k)
            for (j = 0; j < dim.y(); ++j)
              for (i = 0; i < dim.x(); ++i, ++idx1d)
                if (phase_arr[idx1d] == presets::microporous)
                  dpl::sfor<3>([&](auto d) {
                    if (ijk[d] < dim[d] - 1) {
                      auto adj_idx = idx1d + map_idx[d];

                      if (phase_arr[adj_idx] == presets::microporous) {
                        auto new_throat_idx = pnm.throat_count_ + new_throat_incr++;
                      
                        pnm.throat_[attribs::adj][new_throat_idx] = //{0, 0};
                          {old_node_count + voxel_to_row_inc_map[idx1d], old_node_count + voxel_to_row_inc_map[adj_idx]};
                        pnm.throat_[attribs::r_ins][new_throat_idx] = cell_size.x()/4; //min_r_ins;

                        pnm.throat_[attribs::length0][new_throat_idx] = 0;
                        pnm.throat_[attribs::length][new_throat_idx] = cell_size.x();
                        pnm.throat_[attribs::length1][new_throat_idx] = 0;
                      }
                    }
                  });
        }

        std::cout << "\n\nTHROAT_VALUES";


        for (pnm_idx k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j) {
              
            if (auto idx1d = map_idx.dot({0, j, k});
              phase_arr[idx1d] == presets::microporous) {
              auto new_throat_idx = pnm.throat_count_ + new_throat_incr++;
                
              pnm.throat_[attribs::adj][new_throat_idx] = {
                old_node_count + voxel_to_row_inc_map[idx1d],
                pnm.inlet() + darcy_outer_count + darcy_inner_count};
              pnm.throat_[attribs::r_ins][new_throat_idx] = cell_size.x()/4; //min_r_ins;
                
              pnm.throat_[attribs::length0][new_throat_idx] = 0;
              pnm.throat_[attribs::length][new_throat_idx] = cell_size.x();
              pnm.throat_[attribs::length1][new_throat_idx] = 0;
            }

            if (auto idx1d = map_idx.dot({dim.x() - 1, j, k});
              phase_arr[idx1d] == presets::microporous) {
              auto new_throat_idx = pnm.throat_count_ + new_throat_incr++;
                
              pnm.throat_[attribs::adj][new_throat_idx] = {
                old_node_count + voxel_to_row_inc_map[idx1d],
                pnm.outlet() + darcy_outer_count + darcy_inner_count};
              pnm.throat_[attribs::r_ins][new_throat_idx] = cell_size.x()/4; //min_r_ins;
                
              pnm.throat_[attribs::length0][new_throat_idx] = 0;
              pnm.throat_[attribs::length][new_throat_idx] = cell_size.x();
              pnm.throat_[attribs::length1][new_throat_idx] = 0;
            }
          }


          
        pnm.node_count_ += darcy_outer_count + darcy_inner_count;
        pnm.throat_count_ += darcy_outer_count + darcy_darcy_throats + darcy_inlet_outlet;
      }


      if (use_cache && std::filesystem::exists(cache_path)) {
        std::cout << "\n\nUsing CACHED pressure";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(double);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(nrows);
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(double));


        // pnm_idx voxel_to_row_inc = 0;
        //   
        // for (pnm_idx idx1d = 0, k = 0; k < dim.z(); ++k)
        //   for (pnm_idx j = 0; j < dim.y(); ++j)
        //     for (pnm_idx i = 0; i < dim.x(); ++i, ++idx1d)
        //       if (phase_arr[idx1d] == presets::microporous)
        //         voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++;
      }
      else {
        // for (pnm_idx i = 0; i < pnm.node_count_; ++i)
        //   pressure[i] = (1.*i)/pnm.node_count_;

        std::cout << "\n\nparallel_partitioning_START";

        auto partitioning = pnm.parallel_partitioning(processors);
        
        for (auto i = 0; i < processors.prod(); ++i)
          std::cout << std::format("\nblock {}, rows {}--{}, size {}",
            i, partitioning.rows_per_block[i].first, partitioning.rows_per_block[i].second,
            partitioning.rows_per_block[i].second - partitioning.rows_per_block[i].first + 1);

        std::cout << "\n\nparallel_partitioning_END";

        std::cout << "\n\nGeneratePressureInput START...";

        auto input = pnm.generate_pressure_input(partitioning, coef_map);

        std::cout << "\n\nGeneratePressureInput END|||";

        std::cout << "\n\nSave hypre input START...";


        {
          using namespace boost::interprocess;
        
          {
            shared_memory_object smo{open_or_create, "xpm-hypre-input", read_write};
            dpl::hypre::save(input.get_ref(), input.nvalues, partitioning.rows_per_block, smo);
            // input.save(smo);
          }
        
          std::cout << "\n\nSave hypre input END|||";

          std::cout << "\n\nHypre solve MPI START...";
        
          auto start = std::chrono::high_resolution_clock::now();
        


          // {
          //   dpl::hypre::ls_known_ref lk_ref;
          //
          //   shared_memory_object smo{open_only, "xpm-hypre-input", read_only};
          //   mapped_region region(smo, read_only);
          //
          //
          //   std::pair<HYPRE_BigInt, HYPRE_BigInt>* range_ptr;
          //
          //   dpl::hypre::load(region, lk_ref, range_ptr);
          //
          //
          //   auto jlower = 0;
          //   HYPRE_BigInt count = pnm.node_count_;
          //
          //   auto indices = std::make_unique<HYPRE_BigInt[]>(count);
          //   for (HYPRE_BigInt i = 0; i < count; ++i)
          //     indices[i] = jlower + i;
          //   auto pressure_part = std::make_unique<HYPRE_Complex[]>(count);
          //
          //   dpl::hypre::mpi_block::range = {0, pnm.node_count_ - 1};
          //
          //   dpl::hypre::mpi_block::range = *range_ptr;
          //
          //   dpl::hypre::solve(
          //     lk_ref/*input.get_ref()*/,
          //     dpl::hypre::ls_unknown_ref{
          //       count,
          //       indices.get(),
          //       pressure_part.get()
          //     });
          // }



          auto solve_result = std::system(
            std::format("mpiexec -n {} \"{}\" -s",
              processors.prod(),
              "xpm_project.exe"
            ).c_str());


        
          auto stop = std::chrono::high_resolution_clock::now();
        
          cout << "\n\nHypre solve MPI: " <<
            duration_cast<std::chrono::seconds>(stop - start).count() << "s END|||" << endl;
        
        
          shared_memory_object smo{open_only, "xpm-hypre-output", read_only};
          mapped_region mr{smo, read_only};
        
          auto opt_pressure = std::make_unique<double[]>(input.nrows);
        
          std::memcpy(opt_pressure.get(), mr.get_address(), input.nrows*sizeof(double));

          pressure.resize(input.nrows);

          for (HYPRE_BigInt i = 0; i < input.nrows; ++i)
            pressure[partitioning.optimised_to_normal[i]] = opt_pressure[i];

          // auto val = std::accumulate(ptr, ptr + input.nrows, 0.0);
        }

        std::cout << "\n\nPressure solved";


        std::filesystem::create_directory("cache");
        std::ofstream cache_stream(cache_path, std::ofstream::binary);
        cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(double)*pressure.size());
        std::cout << "\n\nPressure cached";
      }


        

      std::vector<pnm_idx> inner_parent(pnm.node_count_);
      std::iota(inner_parent.begin(), inner_parent.end(), 0);
      std::vector<std::uint16_t> inner_rank(pnm.node_count_);
      boost::disjoint_sets inner_ds{inner_rank.data(), inner_parent.data()};
      std::vector<bool> inlet_connected(pnm.node_count_);
      std::vector<bool> outlet_connected(pnm.node_count_);

      auto connected = [&](pnm_idx i) {
        auto rep_set = inner_ds.find_set(i);
        return inlet_connected[rep_set] && outlet_connected[rep_set];
      };

      

      {
        for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
          if (pnm.inner_node(n0) && pnm.inner_node(n1))
            inner_ds.union_set(n0, n1);


        for (auto& [n0, n1] : pnm.throat_.range(attribs::adj)) {
          if (n0 == pnm.inlet())
            inlet_connected[inner_ds.find_set(n1)] = true;
          else if (n1 == pnm.inlet())
            inlet_connected[inner_ds.find_set(n0)] = true;
          else if (n0 == pnm.outlet())
            outlet_connected[inner_ds.find_set(n1)] = true;
          else if (n1 == pnm.outlet())
            outlet_connected[inner_ds.find_set(n0)] = true;
        }

        pnm_idx disconnected_macro = 0;
        for (pnm_idx i = 0; i < old_node_count; ++i)
          if (!connected(i))
            ++disconnected_macro;


        pnm_idx disconnected_darcy = 0;
        for (pnm_idx i = old_node_count; i < pnm.node_count_; ++i)
          if (!connected(i))
            ++disconnected_darcy;


        std::cout << std::format("\n\n Disconnected {} macro, {} darcy nodes", disconnected_macro, disconnected_darcy);
      }



      // std::vector<pnm_idx> inlet_parent;
      // std::vector<std::uint16_t> inlet_rank;
      // std::unique_ptr<boost::disjoint_sets<std::uint16_t*, pnm_idx*>> inlet_ds;
      //
      // std::vector<pnm_idx> outlet_parent;
      // std::vector<std::uint16_t> outlet_rank;
      // std::unique_ptr<boost::disjoint_sets<std::uint16_t*, pnm_idx*>> outlet_ds;
      //
      //
      // std::vector<pnm_idx> total_parent(pnm.node_count_ + 2);
      // std::vector<std::uint16_t> total_rank(pnm.node_count_ + 2);
      // std::iota(total_parent.begin(), total_parent.end(), 0);
      // boost::disjoint_sets total_ds{total_rank.data(), total_parent.data()};

      // {
      //   std::vector<pnm_idx> parent(pnm.node_count_ + 2);
      //   std::vector<std::uint16_t> rank(pnm.node_count_ + 2);
      //
      //   std::iota(parent.begin(), parent.end(), 0);
      //
      //   {
      //     boost::disjoint_sets ds{rank.data(), parent.data()};
      //
      //     for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
      //       if (pnm.inner_node(n0) && pnm.inner_node(n1))
      //         ds.union_set(n0, n1);
      //   }
      //
      //   auto parent_copy = parent;
      //   auto rank_copy = rank;
      //
      //   {
      //     boost::disjoint_sets ds{rank.data(), parent.data()};
      //
      //     for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
      //       if (n0 == pnm.inlet() || n1 == pnm.inlet())
      //         ds.union_set(n0, n1);
      //   }
      //
      //   inlet_parent = std::move(parent);
      //   inlet_rank = std::move(rank);
      //
      //   parent = std::move(parent_copy);
      //   rank = std::move(rank_copy);
      //
      //   {
      //     boost::disjoint_sets ds{rank.data(), parent.data()};
      //
      //     for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
      //       if (n0 == pnm.outlet() || n1 == pnm.outlet())
      //         ds.union_set(n0, n1);
      //   }
      //
      //   outlet_parent = std::move(parent);
      //   outlet_rank = std::move(rank);
      //
      //     
      //
      //   {
      //     // std::vector<pnm_idx> total_parent(pnm.node_count_ + 2);
      //     // std::vector<std::uint16_t> total_rank(pnm.node_count_ + 2);
      //     // std::iota(total_parent.begin(), total_parent.end(), 0);
      //     // boost::disjoint_sets ds{total_rank.data(), total_parent.data()};
      //
      //     for (auto& [n0, n1] : pnm.throat_.range(attribs::adj))
      //       total_ds.union_set(n0, n1);
      //
      //     auto outlet_set = total_ds.find_set(pnm.inlet());
      //
      //     pnm_idx disconnected_macro = 0;
      //     for (pnm_idx i = 0; i < old_node_count; ++i) {
      //       if (total_ds.find_set(i) != outlet_set)
      //         ++disconnected_macro;
      //     }
      //
      //
      //     pnm_idx disconnected_darcy = 0;
      //     for (pnm_idx i = old_node_count; i < pnm.node_count_; ++i)
      //       if (total_ds.find_set(i) != outlet_set)
      //         ++disconnected_darcy;
      //
      //     std::cout << std::format("\n\n Disconnected {} macro, {} darcy nodes [INCORRECT]", disconnected_macro, disconnected_darcy);
      //   }
      //
      //
      //
      //   inlet_ds = std::make_unique<boost::disjoint_sets<std::uint16_t*, pnm_idx*>>(inlet_rank.data(), inlet_parent.data());
      //   outlet_ds = std::make_unique<boost::disjoint_sets<std::uint16_t*, pnm_idx*>>(outlet_rank.data(), outlet_parent.data());
      // }




      {
        double inlet_flow_sum = 0;
        double outlet_flow_sum = 0;

        for (pnm_idx k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j) {
            
            if (auto idx1d = map_idx.dot({0, j, k});
              phase_arr[idx1d] == presets::microporous) {
              if (connected(old_node_count + voxel_to_row_inc_map[idx1d]))
                inlet_flow_sum += -2*cell_size.x()*const_permeability*(1 - pressure[old_node_count + voxel_to_row_inc_map[idx1d]]);
              // else {
              //   // pressure[old_node_count + voxel_to_row_inc_map[idx1d]] = 1;
              //
              //   std::cout << "\nNOT CONNECTED DARCY INLET ADJ TO OUTLET 1";
              // }
            }

            if (auto idx1d = map_idx.dot({dim.x() - 1, j, k});
              phase_arr[idx1d] == presets::microporous) {
              if (connected(old_node_count + voxel_to_row_inc_map[idx1d]))
                outlet_flow_sum += -2*cell_size.x()*const_permeability*(pressure[old_node_count + voxel_to_row_inc_map[idx1d]]);
              // else
              //   std::cout << "\nNOT CONNECTED DARCY OUTLET ADJ TO INLET 2";
            }
          }

        for (pnm_idx i = 0; i < old_throat_count; ++i) {
          auto [n0, n1] = pnm.throat_[attribs::adj][i];

          if (n0 == pnm.inlet()) {
            if (connected(n1))
              inlet_flow_sum += coef_map(i)*(1 - pressure[n1]);
            // else
            //   std::cout << "\nNOT CONNECTED MACRO INLET ADJ TO OUTLET 3";
          }
          else if (n1 == pnm.inlet()) {
            if (connected(n0))
              inlet_flow_sum += coef_map(i)*(1 - pressure[n0]);
            // else
            //   std::cout << "\nNOT CONNECTED MACRO INLET ADJ TO OUTLET 4";
          }
          else if (n0 == pnm.outlet()) {
            if (connected(n1))
              outlet_flow_sum += coef_map(i)*(pressure[n1]);
            // else
            //   std::cout << "\nNOT CONNECTED MACRO OUTLET ADJ TO INLET 5";
          }
          else if (n1 == pnm.outlet()) {
            if (connected(n0))
              outlet_flow_sum += coef_map(i)*(pressure[n0]);
            // else
            //   std::cout << "\nNOT CONNECTED MACRO OUTLET ADJ TO INLET 6";
          }
        }




        
        
        std::cout << std::format("\n\nMICROPOROUS_PERM={} mD\nINLET_PERM={} mD\nOUTLET_PERM={} mD\n",
          const_permeability/darcy_term*1000,
          -inlet_flow_sum/pnm.physical_size.x()/darcy_term*1000,
          -outlet_flow_sum/pnm.physical_size.x()/darcy_term*1000);










        for (auto& [n0, n1] : pnm.throat_.range(attribs::adj)) {
          if (n0 == pnm.inlet() || n0 == pnm.outlet())
            n0 -= pnm.node_count_ - old_node_count;

          if (n1 == pnm.inlet() || n1 == pnm.outlet())
            n1 -= pnm.node_count_ - old_node_count;
        }

        pnm.node_.resize(old_node_count);
        pnm.node_count_ = old_node_count;

        pnm.throat_.resize(old_node_count);
        pnm.throat_count_ = old_throat_count;
      }


      // std::ranges::fill(pressure, 0);

      

      

      









      
      

      

      

      

      

      auto get_pressure = [&pressure](pnm_idx i){ return i < pressure.size() ? pressure[i] : 1; };


      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      assembly->AddPart(CreateNodeActor(pnm, lut_pressure_, get_pressure));
      assembly->AddPart(CreateThroatActor(pnm, lut_pressure_, [&](pnm_idx i) {
        auto [n0, n1] = pnm.throat_[attribs::adj][i];
      
        return (
          (n0 == pnm.inlet() ? 1.0 : n0 == pnm.outlet() ? 0.0 : get_pressure(n0)) +
          (n1 == pnm.inlet() ? 1.0 : n1 == pnm.outlet() ? 0.0 : get_pressure(n1)))/2.0;
      }));

      renderer_->AddActor(assembly);
        
        
        
      



      
      std::cout << "\n\nNetwork actor created";


      


      {
        
          
        {
          auto scale_factor = /*1.0*/pnm.physical_size.x()/dim.x(); // needed for vtk 8.2 floating point arithmetics
            
          img_mapper.Init(scale_factor);
          {
            // auto outlet_set = total_ds.find_set(total_parent.size() - 1);
            auto count = dim.prod<pnm_idx>();

            std::vector<bool> filter(count);
            for (pnm_idx i = 0; i < count; ++i) {

              // auto rep_set = inner_ds.find_set(old_node_count + voxel_to_row_inc_map[i]);


              filter[i] = phase_arr[i] == presets::microporous
               // && inlet_connected[rep_set] && outlet_connected[rep_set]

                // && total_ds.find_set(old_node_count + voxel_to_row_inc_map[i]) == outlet_set
              ;
            }

            img_mapper.Populate(dim, pnm.physical_size/dim, 
              [&](pnm_idx idx1d) { return filter[idx1d]; }
            );
              
            dpl::sfor<6>([&](auto face_idx) {
              GlyphMapperFace& face = std::get<face_idx>(img_mapper.faces_);

              pnm_idx i = 0;
              for (auto idx1d : face.GetIndices()) {
                face.GetColorArray()->SetTypedComponent(i++, 0, 
                  // img_darcy_adj_arr[idx1d]
                  // velem_arr[idx1d].value

                  // total_ds.find_set(old_node_count + voxel_to_row_inc_map[idx1d]) == total_ds.find_set(total_parent.size() - 1) ? 0.5 : 1
                  // 1
                  phase_arr[idx1d] == presets::microporous ? pressure[old_node_count + voxel_to_row_inc_map[idx1d]] : 0

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
              
              // face.POLY_mapper_->SetInputData(face.POLY_polydata__);
              // face.POLY_actor__->SetMapper(face.POLY_mapper_);

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
        0., pnm.physical_size.x(),
        0., pnm.physical_size.y(),
        0., pnm.physical_size.z()};
      tidy_axes_.Build(bounds);

      // tidy_axes_.Build();
      
      // renderer_->ResetCamera();
    }
  };
}