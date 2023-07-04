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
      dpl::sfor<6>([=](auto i) {
        std::get<i>(faces_).Init(half_length);
      });
    }

    void Populate(const v3i& dims, const v3d& cell_size, const auto& filter) {
      auto start = std::chrono::high_resolution_clock::now();
      
      dpl::psfor<6>([=](auto i) {
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

    vtkNew<vtkLookupTable> lut_image_phase_;
    vtkNew<vtkLookupTable> lut_network_element_;



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
    


    bool save_cache = true;

    

    


    // static auto CreateNetworkAssembly(const pore_network_model& pnm, vtkLookupTable* lut) {
    //   auto net = vtkSmartPointer<vtkAssembly>::New();
    //   net->AddPart(CreateNodeActor(pnm, lut, [&](pnm_idx i) { return pnm.node_[attribs::r_ins][i]; }));
    //   net->AddPart(CreateThroatActor(pnm, lut, [&](pnm_idx i) { return pnm.throat_[attribs::r_ins][i]; }));
    //   return net;
    // }




    


    
    
    

    


     
    
  public:
    void LoadImage() {
      auto image_path = R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)";
      auto velems_path = R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)";
      constexpr struct
      {
        std::uint8_t solid = 1;       // dummy value, no '1' is in the image
        std::uint8_t pore = 255;
        std::uint8_t microporous = 0; // we read actual solid '0' as microporous
      } input_spec;


      // auto image_path = R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)";
      // auto velems_path = R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)";
      // constexpr struct
      // {
      //   std::uint8_t solid = 3;       
      //   std::uint8_t pore = 0;
      //   std::uint8_t microporous = 2;  
      // } input_spec;







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
      pnm_3idx dim;

      // std::streamsize

      {
        std::ifstream is(image_path);
        is.seekg(0, std::ios_base::end);

        auto k = is.tellg();
        image_size = is.tellg();
        is.seekg(0, std::ios_base::beg);
        phase_arr = std::make_unique<voxel_tag::phase[]>(image_size);
        is.read(reinterpret_cast<char*>(phase_arr.get()), image_size);

        size_t pore_voxels = 0;
        size_t solid_voxels = 0;
        size_t microporous_voxels = 0;

        for (size_t i = 0; i < image_size; ++i) {
          auto& val = phase_arr[i];
          if (val.value == input_spec.pore) {
            val = presets::pore;
            ++pore_voxels;
          }
          else if (val.value == input_spec.solid) {
            val = presets::solid;
            ++solid_voxels;
          }
          else if (val.value == input_spec.microporous) {
            val = presets::microporous;
            ++microporous_voxels;
          }
        }

        std::cout << std::format("\n\ntotal: {}; pore: {}; solid: {}; microporous: {} voxels", image_size, pore_voxels, solid_voxels, microporous_voxels);
      }

      std::cout << "\n\nImage phases read and array created";

      dim = std::round(std::cbrt(image_size));

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
      auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)";
      v3i dim = 252;

      // auto pnm_path = R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)";
      // v3i dim = 500;

      v3i processors{1};

      auto proc_count = std::thread::hardware_concurrency();

      if (proc_count == 12)
        processors = {2, 2, 3};
      else if (proc_count == 32)
        processors = {4, 4, 2};


      using std::filesystem::path;


      static constexpr auto darcy_term = 9.869233e-13;
      static constexpr auto millidarcy = 10;
      static constexpr auto const_permeability = millidarcy*0.001*darcy_term;

      auto cache_path = std::format("cache/{}-pressure-{:.2f}mD.bin",
        path(pnm_path).parent_path().filename().string(), const_permeability/darcy_term*1e3);



      // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\10x10x10\10x10x10)"
      // R"(C:\dev\.temp\images\SS-1000\XNet)"
      // R"(E:\hwu\research126\d\modelling\networks\TwoScaleNet\MulNet)"

      // auto filter = [this](pnm_idx idx) { return phase_arr[idx] == presets::microporous; };

      
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


      
      lut_network_element_->IndexedLookupOn();
      lut_network_element_->SetNumberOfTableValues(2);
      lut_network_element_->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut_network_element_->SetAnnotation(vtkVariant(0), "Node");
      lut_network_element_->SetTableValue(1, 0.784314, 0.752624, 0.467419);
      lut_network_element_->SetAnnotation(vtkVariant(1), "Throat");


      lut_image_phase_->IndexedLookupOn();  
      lut_image_phase_->SetNumberOfTableValues(3);
      lut_image_phase_->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut_image_phase_->SetAnnotation(vtkVariant(0), "Pore");
      lut_image_phase_->SetTableValue(1, 0.5, 0.5, 0.5);
      lut_image_phase_->SetAnnotation(vtkVariant(1), "Solid");
      lut_image_phase_->SetTableValue(2, 0.8, 0.8, 0.8);
      lut_image_phase_->SetAnnotation(vtkVariant(2), "Microporous");


      // setCentralWidget(qvtk_widget_);


      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      

      pore_network_model pnm{pnm_path, pore_network_model::file_format::statoil};

      

      std::cout << "\n\nNetwork loaded";



      LoadImage();
            
      std::cout << "\n\nLoaded image";

      

      lut_pressure_->SetTableRange(0, 1);


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
      

      if (std::filesystem::exists(cache_path)) {
        std::cout << "\n\nUsing CACHED pressure";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(double);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(nrows);
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(double));


        pnm_idx voxel_to_row_inc = 0;
          
        for (pnm_idx idx1d = 0, k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j)
            for (pnm_idx i = 0; i < dim.x(); ++i, ++idx1d)
              if (phase_arr[idx1d] == presets::microporous)
                voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++;
      }
      else {
        {
          auto image_size = dim.prod<pnm_idx>();

          pnm_idx darcy_outer_count = 0; // darcy nodes connected to macroscopic pore nodes
          pnm_idx darcy_inner_count = 0; // darcy nodes not connected to marcoscopic pore nodes

          for (pnm_idx i = 0; i < image_size; ++i)
            if (phase_arr[i] == presets::microporous)
              if (adj_pore_of_darcy_arr[i] >= 0)
                ++darcy_outer_count;
              else
                ++darcy_inner_count;


          pnm_idx darcy_darcy_throats = 0;
          
          for (pnm_idx k = 1; k < dim.z(); ++k)
            for (pnm_idx j = 1; j < dim.y(); ++j)
              for (pnm_idx i = 1; i < dim.x(); ++i) {
                auto idx1d = map_idx.dot({i, j, k});

                if (phase_arr[idx1d] == presets::microporous) 
                  dpl::sfor<3>([&](auto e) {
                    auto adj_idx = idx1d - map_idx[e];

                    if (phase_arr[adj_idx] == presets::microporous)
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
                    auto new_node_idx = pnm.node_count_ + (voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++);
                    
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
                    auto new_node_idx = pnm.node_count_ + (voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++);
                  
                    pnm.node_[attribs::r_ins][new_node_idx] = cell_size.x()/2;
                    pnm.node_[attribs::pos][new_node_idx] = cell_size*(v3d{i, j, k} + 0.5);
                  }
                }
              }


          std::cout << "\n\nNODES_VALUES";

          for (pnm_idx k = 1; k < dim.z(); ++k)
            for (pnm_idx j = 1; j < dim.y(); ++j)
              for (pnm_idx i = 1; i < dim.x(); ++i) {
                auto idx1d = map_idx.dot({i, j, k});
                
                if (phase_arr[idx1d] == presets::microporous)
                  dpl::sfor<3>([&](auto e) {
                    auto adj_idx = idx1d - map_idx[e];
                    
                    if (phase_arr[adj_idx] == presets::microporous) {
                      auto new_throat_idx = pnm.throat_count_ + new_throat_incr++;
                      
                      pnm.throat_[attribs::adj][new_throat_idx] = //{0, 0};
                        {pnm.node_count_ + voxel_to_row_inc_map[idx1d], pnm.node_count_ + voxel_to_row_inc_map[adj_idx]};
                      pnm.throat_[attribs::r_ins][new_throat_idx] = cell_size.x()/4; //min_r_ins;

                      pnm.throat_[attribs::length0][new_throat_idx] = 0;
                      pnm.throat_[attribs::length][new_throat_idx] = cell_size.x();
                      pnm.throat_[attribs::length1][new_throat_idx] = 0;
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
                  pnm.node_count_ + voxel_to_row_inc_map[idx1d],
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
                  pnm.node_count_ + voxel_to_row_inc_map[idx1d],
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


        std::cout << "\n\nPressure solved";


        std::filesystem::create_directory("cache");
        std::ofstream cache_stream(cache_path, std::ofstream::binary);
        cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(double)*pressure.size());
        std::cout << "\n\nPressure cached";
      }




      
      






      





      

      

      {
        double inlet_flow_sum = 0;
        double outlet_flow_sum = 0;

        for (pnm_idx k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j) {
            
            if (auto idx1d = map_idx.dot({0, j, k});
              phase_arr[idx1d] == presets::microporous) {
              inlet_flow_sum += -2*cell_size.x()*const_permeability*(1 - pressure[pnm.node_count_ + voxel_to_row_inc_map[idx1d]]);
            }

            if (auto idx1d = map_idx.dot({dim.x() - 1, j, k});
              phase_arr[idx1d] == presets::microporous) {
              outlet_flow_sum += -2*cell_size.x()*const_permeability*(pressure[pnm.node_count_ + voxel_to_row_inc_map[idx1d]]);
            }
          }

        for (pnm_idx i = 0; i < pnm.throat_count_; ++i) {
          auto [n0, n1] = pnm.throat_[attribs::adj][i];

          if (n0 == pnm.inlet())
            inlet_flow_sum += coef_map(i)*(1 - pressure[n1]);
          else if (n1 == pnm.inlet())
            inlet_flow_sum += coef_map(i)*(1 - pressure[n0]);
          else if (n0 == pnm.outlet())
            outlet_flow_sum += coef_map(i)*(pressure[n1]);
          else if (n1 == pnm.outlet())
            outlet_flow_sum += coef_map(i)*(pressure[n0]);
         
        }




        
        
        std::cout << std::format("\n\nMICROPOROUS_PERM={} mD\nINLET_PERM={} mD\nOUTLET_PERM={} mD\n",
          const_permeability/darcy_term*1000,
          -inlet_flow_sum/pnm.physical_size.x()/darcy_term*1000,
          -outlet_flow_sum/pnm.physical_size.x()/darcy_term*1000);
      }

      









      
      

      

      

      

      

      auto get_pressure = [&pressure](pnm_idx i){ return i < pressure.size() ? pressure[i] : 1; };


      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      assembly->AddPart(CreateNodeActor(pnm, lut_pressure_, get_pressure));
      assembly->AddPart(CreateThroatActor(pnm, lut_pressure_, [&](pnm_idx i) {
        auto [n0, n1] = pnm.throat_[attribs::adj][i];
      
        return (
          (n0 == pnm.inlet() ? 1.0 : n0 == pnm.outlet() ? 0.0 : get_pressure(n0)) +
          (n1 == pnm.inlet() ? 1.0 : n1 == pnm.outlet() ? 0.0 : get_pressure(n1)))/2.0;
      }));


      // assembly->AddPart(CreateNodeActor(pnm, lut_network_element_, [](pnm_idx i) { return 0; }));
      // assembly->AddPart(CreateThroatActor(pnm, lut_network_element_, [](pnm_idx i) { return 1; }));

      renderer_->AddActor(assembly);
        
        
        
      



      
      std::cout << "\n\nNetwork actor created";


      


      {
        
          
        {
          auto scale_factor = /*1.0*/pnm.physical_size.x()/dim.x(); // needed for vtk 8.2 floating point arithmetics
            
          img_mapper.Init(scale_factor);
          {
            img_mapper.Populate(dim, pnm.physical_size/dim, 
              [this](pnm_idx idx1d) {
                // return true;
                return phase_arr[idx1d] == presets::microporous;
              }
            );
              
            dpl::sfor<6>([&](auto face_idx) {
              GlyphMapperFace& face = std::get<face_idx>(img_mapper.faces_);

              pnm_idx i = 0;
              for (auto idx1d : face.GetIndices()) {
                face.GetColorArray()->SetTypedComponent(i++, 0, 
                  // img_darcy_adj_arr[idx1d]
                  // velem_arr[idx1d].value

                  phase_arr[idx1d] == presets::microporous ? pressure[pnm.node_count_ + voxel_to_row_inc_map[idx1d]] : 0

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