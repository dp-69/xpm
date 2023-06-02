#pragma once


#include "xpm/functions.h"

#include <dpl/hypre/Input.hpp>
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
  namespace geometric_properties
  {
    struct equilateral_triangle_properties
    {
      static constexpr double area(double r_ins = 1) {
        return 5.19615242271*r_ins*r_ins;
      }

      // k * G, k - coefficient, G - shape factor
      static constexpr double conductance(double area = 1, double viscosity = 1) {
        return 0.0288675134595*area*area/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }
    };
  }


  
  template<int face_idx>
  class ImageDataGlyphMapperFace
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

    vtkNew<vtkPoints> points_;
    
  public:
    vtkNew<vtkGlyph3DMapper> glyphs_;
    vtkNew<vtkIntArray> darcy_adj_out_;
    vtkNew<vtkActor> actor_;
    
    void Init(double half_length = 0.5) {
      darcy_adj_out_->SetName("darcy_adj");

      vtkNew<vtkPolyData> polydata;
      polydata->GetPointData()->SetScalars(darcy_adj_out_);
      polydata->SetPoints(points_);

      glyphs_->OrientOff();
      glyphs_->SetScaleFactor(half_length); 
      glyphs_->SetScaleModeToNoDataScaling();
      glyphs_->SetInputData(polydata);
      glyphs_->SetSourceData(Quad());
    }

    void Populate(const v3i& cells, const v3d& cell_size, const auto& filter, const auto& post) {
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
              post(idx1d);
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

              if constexpr (face::is_upper) {
                if (filtered) {
                  points_->InsertNextPoint(pos);
                  post(idx1d);
                }
              }
              else {
                if (adj_filtered) {
                  points_->InsertNextPoint(pos);
                  post(adj_idx1d);
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
              post(idx1d);
            }
          }
        }
    }
  };


  class ImageDataGlyphMapper
  {
  public:

    void Init(double half_length) {
      dpl::sfor<6>([this, half_length](auto i) {
        std::get<i>(faces_).Init(half_length);
      });
    }

    template<typename Filter>
    void Populate(const v3i& dims, const v3d& cell_size, const Filter& filter, vtkIntArray* velems_adj_arr) {
      auto start = std::chrono::high_resolution_clock::now();
      
      dpl::psfor<6>([=, this](auto i) {
        auto* face_mapper = &std::get<i>(faces_);

        face_mapper->Populate(dims, cell_size, filter, 
          [=](pnm_idx idx) {
            auto val = velems_adj_arr->GetTypedComponent(idx, 0);
            face_mapper->darcy_adj_out_->InsertNextTypedTuple(&val);
          }
        );

        // std::cout << "\n\nFaces " << i;
      });

      auto stop = std::chrono::high_resolution_clock::now();
 
      cout << "\n\nFaces total time: " <<
        duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << endl;
    }

    std::tuple<
      ImageDataGlyphMapperFace<0>,
      ImageDataGlyphMapperFace<1>,
      ImageDataGlyphMapperFace<2>,
      ImageDataGlyphMapperFace<3>,
      ImageDataGlyphMapperFace<4>,
      ImageDataGlyphMapperFace<5>
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
    vtkNew<vtkLookupTable> lut_pore_solid_;
    vtkNew<vtkLookupTable> lut_velem_;

    
    
    
    vtkNew<vtkImageData> image_data_;
    // vtkNew<vtkThreshold> threshold_;



    // vtkNew<vtkActor> image_actor_;




    ImageDataGlyphMapper img_mapper;




    vtkNew<vtkIntArray> velem_array;
    
    
    /*
     * connected pore node of a solid voxel. Gives a cluster number [0, n].
     * -2 for pore voxels,
     * -1 for solid voxels not connected
     * 
     */
    vtkNew<vtkIntArray> img_darcy_adj_array;


    
    

    


    static auto CreateNetworkAssembly(const pore_network_model& pnm, vtkLookupTable* lut) {
      auto net = vtkSmartPointer<vtkAssembly>::New();
      net->AddPart(CreateNodeActor(pnm, lut, [&](pnm_idx i) { return pnm.node_[attribs::r_ins][i]; }));
      net->AddPart(CreateThroatActor(pnm, lut, [&](pnm_idx i) { return pnm.throat_[attribs::r_ins][i]; }));
      return net;
    }




    


    
    
    

    



    
  public:
    void LoadImage() {
      auto image_path = 
        // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
        
        R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)"
        // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\pnm_petronas\images\Est_3phase500cubed4micron_NORM.raw)"
      ;

      auto velems_path = 
        R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
        // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
      ;




      




      

      


      
      lut_pore_solid_->IndexedLookupOn();
      
      lut_pore_solid_->SetNumberOfTableValues(2);
      

      lut_pore_solid_->SetTableValue(0, 0.6, 0.6, 0.6);
      lut_pore_solid_->SetAnnotation(vtkVariant(0), "SOLID");

      lut_pore_solid_->SetTableValue(1, 0, 0, 0);
      lut_pore_solid_->SetAnnotation(vtkVariant(1), "PORE");
      
      // lut_pore_solid_->SetTableRange(0, 230);
      // lut_pore_solid_->Modified();


      vtkNew<vtkUnsignedCharArray> phase_array;
      


      

      
      
      size_t size;
      dpl::vector_n<pnm_idx, 3> dim;


      std::vector<unsigned char> phase_buffer;
      
      {
        // Open the stream
        std::ifstream is(image_path);
        // Determine the file length
        is.seekg(0, std::ios_base::end);
        size = is.tellg();
        is.seekg(0, std::ios_base::beg);
        phase_buffer.resize(size);
        is.read(reinterpret_cast<char*>(phase_buffer.data()), size);


        
        for (size_t i = 0; i < size; ++i)
          phase_array->InsertTypedComponent(i, 0, phase_buffer[i]/*v[i] == 0 ? 0 : 1*/);

        
          // phase_array->InsertTypedComponent(i, 0, v[i] == 0 ? 230 : 99);
          // phase_array->InsertNextTuple1(/*i<(size/2) ? 1 : 0*/v[i] == 0 ? 230 : 99);

        phase_array->SetName("phase");
        // image_data_->GetCellData()->SetScalars(phase_array);
        image_data_->GetCellData()->AddArray(phase_array);
      }

      std::cout << "\n\nImage phases read and array created";

      dim = std::round(std::cbrt(size));

      {
        auto velems = pore_network_model::read_icl_velems(velems_path, dim);

        std::cout << "\n\nVelems file read";

        // auto [minelem, maxelem] = std::ranges::minmax_element(velems);
        
        for (pnm_idx i = 0, count = velems.size(); i < count; ++i) {
          // if (phase_buffer[i] == 3 && velems[i] != -2) {
          //   int error = 1;
          // }
          
          velem_array->InsertTypedComponent(i, 0, velems[i]);
        }    
        
        
        velem_array->SetName("velem");
        image_data_->GetCellData()->AddArray(velem_array);

        std::cout << "\n\nVelems array filled";

        
        auto [min, max] = std::ranges::minmax_element(velems);
        auto count = *max - *min + 1;

        std::cout << "\n\nRANGES MINMAX VELEMS";
        
        // lut_velem_->IndexedLookupOn();
        lut_velem_->SetNumberOfTableValues(count);

        // vtkStdString nan_text = "Nan";
        
        for (int32_t i = 0; i < count; ++i) {
          auto coef = static_cast<double>(i*45%count)/count;

          auto color = QColor::fromHsl(coef*255, 175, 122);
          lut_velem_->SetTableValue(i, color.redF(), color.greenF(), color.blueF());  // NOLINT(clang-diagnostic-double-promotion)
          // lut_velem_->SetAnnotation(vtkVariant(i), nan_text/*std::to_string(i)*/); // KILLS PERFORMANCE A LOT!

          
        }

        // lut_velem_->ResetAnnotations();
        lut_velem_->SetTableValue(0, 0.4, 0.4, 0.4);
        lut_velem_->SetTableValue(1, 0.55, 0.55, 0.55);
        lut_velem_->SetTableRange(*min, *max);
        

        std::cout << "\n\nlut_velem_ created";
        
        // dpl::vtk::PopulateLutRedWhiteBlue(lut_velem_);
        // lut_velem_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);

        {
          img_darcy_adj_array->SetNumberOfComponents(1);
          // img_darcy_adj_array->SetNumberOfTuples(dim.prod());
          
          pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};
        
          pnm_idx idx1d = 0;
          
          for (pnm_idx k = 0; k < dim.z(); ++k)
            for (pnm_idx j = 0; j < dim.y(); ++j)
              for (pnm_idx i = 0; i < dim.x(); ++i, ++idx1d) {
                int32_t adj = -1; // Solid not connected 
                
                if (velems[idx1d] < 0) { // Solid
                  if (i > 0)
                    if (auto adj_velem = velems[idx1d - map_idx.x()]; adj_velem > 1)
                      adj = adj_velem;
        
                  if (i < dim.x() - 1)
                    if (auto adj_velem = velems[idx1d + map_idx.x()]; adj_velem > 1)
                      adj = adj_velem;
        
        
                  if (j > 0)
                    if (auto adj_velem = velems[idx1d - map_idx.y()]; adj_velem > 1)
                      adj = adj_velem;
        
                  if (j < dim.y() - 1)
                    if (auto adj_velem = velems[idx1d + map_idx.y()]; adj_velem > 1)
                      adj = adj_velem;
        
        
                  if (k > 0)
                    if (auto adj_velem = velems[idx1d - map_idx.z()]; adj_velem > 1)
                      adj = adj_velem;
        
                  if (k < dim.z() - 1)
                    if (auto adj_velem = velems[idx1d + map_idx.z()]; adj_velem > 1)
                      adj = adj_velem;
                }
                else
                  adj = -2; // Pore voxel

                img_darcy_adj_array->InsertNextTypedTuple(&adj);
              }
        
          
          img_darcy_adj_array->SetName("darcy_adj");
          image_data_->GetCellData()->AddArray(img_darcy_adj_array);
        
          std::cout << "\n\nVelems_adj array produced and filled";
        }
      }


      

      
      
      
      auto dim_side = std::round(std::cbrt(size));
      image_data_->SetDimensions(v3i{dim_side + 1});
      image_data_->SetSpacing(v3d{1.0/dim_side});
      image_data_->SetOrigin(v3d{0});

      
      


      vtkDataArray* selected_arr = img_darcy_adj_array;

      // vtkDataArray* selected_arr = velem_array;
      
      
      image_data_->GetCellData()->SetActiveScalars(selected_arr->GetName());

      // vtkNew<vtkDataSetMapper> mapper;
      //
      // {
      //   mapper->SetInputData(image_data_);
      //   mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, selected_arr->GetName()); // This is needed only for vtkImageData without vtkThreshold
      //
      //   mapper->SetColorModeToMapScalars(); 
      //   mapper->UseLookupTableScalarRangeOn();
      //   mapper->SetScalarModeToUseCellData();
      //   mapper->SetLookupTable(image_data_->GetCellData()->GetScalars() == phase_array ? lut_pore_solid_ : lut_velem_);
      // }
      
      // {
      //
      //   threshold_->SetInputData(image_data_);
      //   threshold_->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, velem_adjacent_array->GetName());
      //   
      //   // threshold_->SetLowerThreshold(2);
      //   // threshold_->SetUpperThreshold(100000);
      //   
      //   threshold_->SetLowerThreshold(5);
      //   threshold_->SetUpperThreshold(1e9);
      //
      //
      //   mapper->SetInputData(threshold_->GetOutput());
      //   
      //   // mapper->SetInputData(image_data_);
      //   // mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, selected_arr->GetName()); // This is needed only for vtkImageData without vtkThreshold
      //
      //   mapper->SetColorModeToMapScalars(); 
      //   mapper->UseLookupTableScalarRangeOn();
      //   mapper->SetScalarModeToUseCellData();
      //   mapper->SetLookupTable(image_data_->GetCellData()->GetScalars() == phase_array ? lut_pore_solid_ : lut_velem_);
      //
      //   image_data_->Modified();
      //   threshold_->Modified();
      //   threshold_->Update();
      //
      // }
      



      
      
      

      
      //
      // image_actor_->SetMapper(mapper);
      // image_actor_->GetProperty()->SetEdgeVisibility(false/*true*/);
      // image_actor_->GetProperty()->SetEdgeColor(v3d{0.5} /*0, 0, 0*/);
      //
      // image_actor_->GetProperty()->SetAmbient(0.5);
      // image_actor_->GetProperty()->SetDiffuse(0.4);
      // image_actor_->GetProperty()->BackfaceCullingOn();
      

      // renderer_->AddActor(image_actor_);
    }


    auto SolvePressure(pore_network_model& pnm) {

      using eq_tri = geometric_properties::equilateral_triangle_properties;
      using namespace attribs;
      
      // auto calc_coef = [&pnm](pnm_idx i) {
      //   auto [n0, n1] = pnm.throat_[adj][i];
      //   
      //   return
      //     (pnm.inner_node(n0) ? pnm.throat_[length0][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n0])) : 0.0) +
      //     (pnm.inner_node(n1) ? pnm.throat_[length1][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n1])) : 0.0) +
      //     pnm.throat_[length][i]/eq_tri::conductance(eq_tri::area(pnm.throat_[r_ins][i]));
      // };

      dpl::hypre::SparseMatrix matrix(pnm.node_count_);

      std::vector<double> free_terms(pnm.node_count_, 0);
      
      
      for (pnm_idx i = 0; i < pnm.throat_count_; ++i) {
        auto [n0, n1] = pnm.throat_[adj][i];

        auto k = pnm.throat_[length0][i];
        auto w = pnm.throat_[length1][i];
        auto Q = pnm.throat_[length][i];
        
        
        auto coef = -1.0/(
          (pnm.inner_node(n0) ? pnm.throat_[length0][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n0])) : 0.0) +
          (pnm.inner_node(n1) ? pnm.throat_[length1][i]/eq_tri::conductance(eq_tri::area(pnm.node_[r_ins][n1])) : 0.0) +
          pnm.throat_[length][i]/eq_tri::conductance(eq_tri::area(pnm.throat_[r_ins][i])));

        if (n0 == pnm.inlet()) {
          free_terms[n1] += coef/**1 Pa*/;
          matrix.AddDiagCoef(n1, coef);
        }
        else if (n1 == pnm.inlet()) {
          free_terms[n0] += coef/**1 Pa*/;
          matrix.AddDiagCoef(n0, coef);
        }
        else if (n0 == pnm.outlet()) {
          // free_terms[n1] += coef/**0 Pa*/;
          matrix.AddDiagCoef(n1, coef);
        }
        else if (n1 == pnm.outlet()) {
          // free_terms[n0] += coef/**0 Pa*/;
          matrix.AddDiagCoef(n0, coef);
        }
        else /*if (n0 != pnm.outlet() && n1 != pnm.outlet())*/ {
          matrix.AddDifferenceCoefs(n0, n1, -coef);
        }
      }


      dpl::hypre::Input input{matrix, std::move(free_terms)};

      auto values = input.Solve();

      return values;

      // int p = 3;
      
    }
    
    void Init() {
      auto pnm_path = 
        // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\10x10x10\10x10x10)"
        // R"(C:\dev\.temp\images\SS-1000\XNet)"
        // R"(E:\hwu\research126\d\modelling\networks\TwoScaleNet\MulNet)"

        R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
        // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
      ;

      v3i dim = 252;
      // v3i dim = 500;

      vtkUnsignedCharArray* phase_in;
      // vtkIntArray* velems_arr_in;
      // vtkIntArray* darcy_adj_in;

      auto filter = [&, this](pnm_idx idx) {
        // return true;

        // vtkIntArray* velems_adj_arr_in;

        // return velems_adj_arr_in->GetTypedComponent(idx, 0) > 2;

        
        // return phase_in->GetTypedComponent(idx, 0) == 2;

        return velem_array->GetTypedComponent(idx, 0) == -2;

        // ;/

        int z = idx / (dim.x() * dim.y());
        int y = (idx - z * dim.x() * dim.y()) / dim.x();
        int x = idx - z * dim.x() * dim.y() - y * dim.x();

        // int z = i%dim.z();
        // int y = (i/dim.z())%dim.y();
        // int x = i/(dim.y()*dim.z()); 
      };

     
      

      
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
          return std::get<0>(img_mapper.faces_).actor_->GetProperty()->GetEdgeVisibility();
        };
        edges.set = [this](bool v) {
          dpl::sfor<6>([this, v](auto i) {
            std::get<i>(img_mapper.faces_).actor_->GetProperty()->SetEdgeVisibility(v);
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


      

     



      


      

      



      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      

      pore_network_model pnm{pnm_path, pore_network_model::file_format::statoil};

      std::cout << "\n\nNetwork loaded";


      


      LoadImage();
            
      std::cout << "\n\nLoaded image";




      












      

     

      lut_pressure_->SetTableRange(0, 1);


      using namespace attribs;
      // static constexpr const auto& pos = attribs::pos;
      // using pos = attribs::pos;




      


      
      {
        // auto pred = [](int32_t darcy_adj) {
        //   // return darcy_adj == 10;
        //   
        //   return/* darcy_adj == 9;*/ darcy_adj >= 0/* && darcy_adj <= 0*/;
        // };

        
        auto min_r_ins = *std::ranges::min_element(pnm.throat_.range(r_ins));

        

        // int max_solid_non = 1000;
        
        
        auto* ptr = static_cast<int*>(img_darcy_adj_array->GetVoidPointer(0));

        auto darcy_node_con_cluster_adj_count =
          std::count_if(ptr, ptr + img_darcy_adj_array->GetNumberOfTuples(), [&](int32_t darcy_adj) { return darcy_adj >= 0; });

        auto darcy_node_inner_adj_count =
          std::count_if(ptr, ptr + img_darcy_adj_array->GetNumberOfTuples(), [&](int32_t darcy_adj) { return darcy_adj == -1; });
        
      
        
        auto* velems = static_cast<int*>(velem_array->GetVoidPointer(0));
        pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};


        pnm_idx darcy_darcy_throats = 0;


        std::unordered_map<pnm_idx, pnm_idx> voxel_to_row_inc_map;
        
        for (pnm_idx k = 1; k < dim.z(); ++k)
          for (pnm_idx j = 1; j < dim.y(); ++j)
            for (pnm_idx i = 1; i < dim.x(); ++i)
              if (auto idx1d = map_idx.dot({i, j, k});
                velems[idx1d] == -2) // Solid
                dpl::sfor<3>([&](auto e) {
                  if (velems[idx1d - map_idx[e]] == -2)
                    ++darcy_darcy_throats;  
                });


        darcy_darcy_throats = 0;

        
        

        // darcy_adj_count = 0;

        pnm.node_.resize(pnm.node_count_ + darcy_node_con_cluster_adj_count + darcy_node_inner_adj_count);
        pnm.throat_.resize(pnm.throat_count_ + darcy_node_con_cluster_adj_count + darcy_darcy_throats);
        
        
        pnm_idx new_throat_incr = 0;
        
        auto cell_size = pnm.physical_size/dim;


        
        for (auto& [n0, n1] : pnm.throat_.range(attribs::adj)) {
          if (n0 == pnm.inlet() || n0 == pnm.outlet())
            n0 += darcy_node_con_cluster_adj_count + darcy_node_inner_adj_count;

          if (n1 == pnm.inlet() || n1 == pnm.outlet())
            n1 += darcy_node_con_cluster_adj_count + darcy_node_inner_adj_count;
        }

        
        pnm_idx voxel_to_row_inc = 0;
        
        for (pnm_idx idx1d = 0, k = 0; k < dim.z(); ++k)
          for (pnm_idx j = 0; j < dim.y(); ++j)
            for (pnm_idx i = 0; i < dim.x(); ++i, ++idx1d) {
              if (ptr[idx1d] >= 0) { // Solid
                voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++;


                auto new_node_idx = pnm.node_count_ + voxel_to_row_inc_map[idx1d];
                
                pnm.node_[attribs::r_ins][new_node_idx] = cell_size.x()/2;
                pnm.node_[attribs::pos][new_node_idx] = cell_size*(v3d{i, j, k} + 0.5);
                

                auto new_throat_idx = pnm.throat_count_ + new_throat_incr;
                
                pnm.throat_[adj][new_throat_idx] = {new_node_idx, ptr[idx1d]};
                pnm.throat_[r_ins][new_throat_idx] = cell_size.x()/4;//min_r_ins;
                
                pnm.throat_[length0][new_throat_idx] = 0;
                pnm.throat_[length][new_throat_idx] =
                  (pnm.node_[attribs::pos][new_node_idx] - pnm.node_[attribs::pos][ptr[idx1d]]).length();
                pnm.throat_[length1][new_throat_idx] = 0;
                
                               
                ++new_throat_incr;
              }
              else if (ptr[idx1d] == -1) {
                voxel_to_row_inc_map[idx1d] = voxel_to_row_inc++;
              
                auto new_node_idx = pnm.node_count_ + voxel_to_row_inc_map[idx1d];
              
                pnm.node_[attribs::r_ins][new_node_idx] = cell_size.x()/2;
                pnm.node_[attribs::pos][new_node_idx] = cell_size*(v3d{i, j, k} + 0.5);
              }
              

              
              
            }


        pnm_idx loaded = 0;
        
        for (pnm_idx k = 1; k < dim.z(); ++k)
          for (pnm_idx j = 1; j < dim.y(); ++j)
            for (pnm_idx i = 1; i < dim.x(); ++i) {
              // int32_t adj = -1; // Solid not connected 

              auto idx1d = map_idx.dot({i, j, k});
              
              if (velems[idx1d] == -2) {
                // Solid

                dpl::sfor<3>([&](auto e) {
                  auto adj_idx = idx1d - map_idx[e];
                  if (loaded < darcy_darcy_throats && velems[adj_idx] == -2) {
                    auto new_throat_idx = pnm.throat_count_ + new_throat_incr;
                    
                    pnm.throat_[adj][new_throat_idx] = //{0, 0};
                      {pnm.node_count_ + voxel_to_row_inc_map[idx1d], pnm.node_count_ + voxel_to_row_inc_map[adj_idx]};
                    pnm.throat_[r_ins][new_throat_idx] = cell_size.x()/4; //min_r_ins;

                    pnm.throat_[length0][new_throat_idx] = 0;
                    pnm.throat_[length][new_throat_idx] = cell_size.x();
                    pnm.throat_[length1][new_throat_idx] = 0;

                    ++new_throat_incr;
                    ++loaded;
                  }
                });
              }
            }


        
        
        pnm.node_count_ += darcy_node_con_cluster_adj_count + darcy_node_inner_adj_count;

        pnm.throat_count_ += darcy_node_con_cluster_adj_count + darcy_darcy_throats;

        
      
      
        
        
        
        
      }


      auto pressure = SolvePressure(pnm);

      std::cout << "\n\nPressure solved";
      


      
      

      

      

      

      auto get_pressure = [&pressure](pnm_idx i){ return i < pressure.size() ? pressure[i] : 1; };
      
      
        auto assembly = vtkSmartPointer<vtkAssembly>::New();
        assembly->AddPart(CreateNodeActor(pnm, lut_pressure_, 
          get_pressure
        //   [&](pnm_idx i) {
        //   return i < pressure.size() ? pressure[i] : 0;
        //
        //   
        // }
        
        ));
        assembly->AddPart(CreateThroatActor(pnm, lut_pressure_, [&](pnm_idx i) {
          auto [n0, n1] = pnm.throat_[adj][i];

          return (
            (n0 == pnm.inlet() ? 1.0 : n0 == pnm.outlet() ? 0.0 : get_pressure(n0)/*pressure[n0]*/) +
            (n1 == pnm.inlet() ? 1.0 : n1 == pnm.outlet() ? 0.0 : get_pressure(n1)/*pressure[n1]*/))/2.0;

          // return 0;pressure[i];
        }));
        renderer_->AddActor(assembly);
        
        
        
        
        // color_array->SetName("color");
        
        
        
        
        
        std::cout << "\n\nNetwork actor created";

      
      
      


      if (false) {
        // auto [min, max] = std::ranges::minmax_element(icl_pnm_inv.node_.range(attribs::r_ins));
        //
        // lut_continuous_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);
        //
        // renderer_->AddActor(CreateNetworkAssembly(icl_pnm_inv, lut_continuous_));
        //
        // std::cout << "\n\nNetwork actor created";


        
        
       
      }
      else
      {

        // auto [min, max] = std::ranges::minmax_element(pnm.node_.range(attribs::r_ins));
        //
        // lut_continuous_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);
        //
        // renderer_->AddActor(CreateNetworkAssembly(pnm, lut_continuous_));
        //
        // std::cout << "\n\nNetwork actor created";
        


        
        

        
        // image_actor_->SetUserTransform()
        // {
        //   vtkNew<vtkTransform> trans;
        //   trans->PostMultiply();
        //   trans->Scale(v3d{pnm.physical_size.x()});
        //   // trans->Translate(icl_pnm_inv.physical_size.x(), 0, 0);
        //   // image_actor_->SetUserTransform(trans);
        // }


        {
          phase_in = static_cast<vtkUnsignedCharArray*>(image_data_->GetCellData()->GetArray("phase"));
          // velems_arr_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("velem"));
          // darcy_adj_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("darcy_adj"));
          
          {
            auto scale_factor = /*1.0*/pnm.physical_size.x()/dim.x(); // needed for vtk 8.2 floating point arithmetics
            
            
            
            img_mapper.Init(scale_factor);
            



            {
              img_mapper.Populate(dim, pnm.physical_size/dim, filter, img_darcy_adj_array);
              
              dpl::sfor<6>([&](auto i) {
                auto& mapper = std::get<i>(img_mapper.faces_);
                
                vtkGlyph3DMapper* glyphs = mapper.glyphs_.Get();
                vtkActor* actor = mapper.actor_.Get();



                glyphs->SetLookupTable(lut_velem_);
                glyphs->SetColorModeToMapScalars();
                glyphs->UseLookupTableScalarRangeOn();
                glyphs->SetScalarModeToUsePointData();


                actor->SetMapper(glyphs);

                actor->GetProperty()->SetEdgeVisibility(/*false*/false);
                actor->GetProperty()->SetEdgeColor(v3d{0.25} /*0, 0, 0*/);
                
                actor->GetProperty()->SetAmbient(0.5);
                actor->GetProperty()->SetDiffuse(0.4);
                actor->GetProperty()->BackfaceCullingOn();
                
                // renderer_->AddActor(actor);
              });
            }
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