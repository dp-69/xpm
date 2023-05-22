#pragma once


#include "xpm/functions.h"

#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>
#include <dpl/qt/property_editor/QPropertyTreeView.hpp>
  
// #include <QWidget>
#include <QMainWindow>
#include <QSplitter>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkCylinderSource.h>
#include <vtkFieldData.h>
#include <vtkCellData.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGlyph3D.h>
#include <vtkGlyph3DMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkConeSource.h>
#include <vtkVersionMacros.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkDataSetMapper.h>
#include <vtkTransform.h>
#include <vtkThreshold.h>
#include <vtkAssembly.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>

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
    
    std::unique_ptr<QPropertyTreeView> tree_view_factory_;
    QPropertyTreeView* tree_view_ = nullptr; 

    
    
    vtkNew<vtkRenderer> renderer_; 
    vtkNew<vtkGenericOpenGLRenderWindow> render_window_;
    
    vtkRenderWindowInteractor* interactor_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;
    vtkNew<vtkLookupTable> lut_continuous_;
    vtkNew<vtkLookupTable> lut_pore_solid_;
    vtkNew<vtkLookupTable> lut_velem_;

    
    
    
    vtkNew<vtkImageData> image_data_;
    // vtkNew<vtkThreshold> threshold_;



    vtkNew<vtkActor> image_actor_;




    


    
    

    


    static auto CreateNetworkAssembly(const pore_network_model& pnm, vtkLookupTable* lut) {
      auto net = vtkSmartPointer<vtkAssembly>::New();
      net->AddPart(CreateNodeActor(pnm, lut));
      net->AddPart(CreateThroatActor(pnm, lut));
      return net;
    }




    template<int e1_idx, typename Filter, typename Post>
    void FilterFaces(const v3i& dims, const v3d& cell_size, const Filter& filter, const Post& post, vtkPolyData* polydata) {
      static constexpr auto e1_dim = dpl::sdim<3, e1_idx>{};
      static constexpr auto e2_dim = e1_dim.next();
      static constexpr auto e3_dim = e2_dim.next();

      pnm_3idx map_idx{1, dims.x(), dims.x()*dims.y()};
      pnm_3idx ijk;

      auto& e3 = ijk[e3_dim];
      auto& e2 = ijk[e2_dim];
      auto& e1 = ijk[e1_dim];

      auto e3_count = dims[e3_dim];
      auto e2_count = dims[e2_dim];
      auto e1_count = dims[e1_dim];
      
      auto adj_step = map_idx[e1_dim];

      v3d pos;

      auto* points = polydata->GetPoints();
      auto* polys = polydata->GetPolys();

      auto point_count = points->GetNumberOfPoints();

      pnm_idx idx1d;
      
      for (e3 = 0; e3 < e3_count; ++e3)
        for (e2 = 0; e2 < e2_count; ++e2) {
          e1 = 0;

          idx1d = ijk.dot(map_idx);
          
          if (filter(idx1d)) {
            pos[e1_dim] = (0)*cell_size[e1_dim];
          
            pos[e2_dim] = (e2)*cell_size[e2_dim];
            pos[e3_dim] = (e3)*cell_size[e3_dim];
            points->InsertNextPoint(pos);
          
            pos[e2_dim] = (e2 + 1)*cell_size[e2_dim];
            points->InsertNextPoint(pos);
          
            pos[e3_dim] = (e3 + 1)*cell_size[e3_dim];
            points->InsertNextPoint(pos);
          
            pos[e2_dim] = (e2)*cell_size[e2_dim];
            points->InsertNextPoint(pos);
          
            polys->InsertNextCell({point_count + 3, point_count + 2, point_count + 1, point_count + 0});
            post(idx1d);
            point_count += 4;
          }
          
          for (; e1 < e1_count - 1; ++e1) {
            idx1d = ijk.dot(map_idx);
            bool idx1d_filter = filter(idx1d);

            auto adj_idx1d = idx1d + adj_step;
            bool adj_idx1d_filter = filter(adj_idx1d);

            if (idx1d_filter != adj_idx1d_filter) {
              pos[e1_dim] = (e1 + 1)*cell_size[e1_dim];

              pos[e2_dim] = (e2)*cell_size[e2_dim];
              pos[e3_dim] = (e3)*cell_size[e3_dim];
              points->InsertNextPoint(pos);

              pos[e2_dim] = (e2 + 1)*cell_size[e2_dim];
              points->InsertNextPoint(pos);

              pos[e3_dim] = (e3 + 1)*cell_size[e3_dim];
              points->InsertNextPoint(pos);

              pos[e2_dim] = (e2)*cell_size[e2_dim];
              points->InsertNextPoint(pos);

              if (idx1d_filter) {
                polys->InsertNextCell({point_count, point_count + 1, point_count + 2, point_count + 3});
                post(idx1d);
              }
              else {
                polys->InsertNextCell({point_count + 3, point_count + 2, point_count + 1, point_count + 0});
                post(adj_idx1d);
              }

              point_count += 4;
            }
          }


          idx1d = ijk.dot(map_idx);
          if (filter(idx1d)) {
            pos[e1_dim] = (e1_count)*cell_size[e1_dim];
          
            pos[e2_dim] = (e2)*cell_size[e2_dim];
            pos[e3_dim] = (e3)*cell_size[e3_dim];
            points->InsertNextPoint(pos);
          
            pos[e2_dim] = (e2 + 1)*cell_size[e2_dim];
            points->InsertNextPoint(pos);
          
            pos[e3_dim] = (e3 + 1)*cell_size[e3_dim];
            points->InsertNextPoint(pos);
          
            pos[e2_dim] = (e2)*cell_size[e2_dim];
            points->InsertNextPoint(pos);
          
            polys->InsertNextCell({point_count, point_count + 1, point_count + 2, point_count + 3});
            post(idx1d);
            point_count += 4;
          }
        }
    }


    template<int e1_idx, typename Filter, typename Post>
    void FilterFacesGlyph(
      const v3i& dims, const v3d& cell_size, const Filter& filter, const Post& post, vtkFloatArray* orient_array, vtkPoints* points) {

      static constexpr auto e1_dim = dpl::sdim<3, e1_idx>{};
      static constexpr auto e2_dim = e1_dim.next();
      static constexpr auto e3_dim = e2_dim.next();

      pnm_3idx map_idx{1, dims.x(), dims.x()*dims.y()};
      pnm_3idx ijk;

      auto& e3 = ijk[e3_dim];
      auto& e2 = ijk[e2_dim];
      auto& e1 = ijk[e1_dim];

      auto e3_count = dims[e3_dim];
      auto e2_count = dims[e2_dim];
      auto e1_count = dims[e1_dim];
      
      auto adj_step = map_idx[e1_dim];

      v3d pos;
      dpl::vector_n<float, 3> orient_vec{0};
      dpl::vector_n<float, 3> orient_vec_opp{0};
      orient_vec[e1_dim] = 1;
      orient_vec_opp[e1_dim] = -1;

      // auto* points = polydata->GetPoints();
      // auto* polys = polydata->GetPolys();

      pnm_idx idx1d;
      
      for (e3 = 0; e3 < e3_count; ++e3)
        for (e2 = 0; e2 < e2_count; ++e2) {
          e1 = 0;

          idx1d = ijk.dot(map_idx);
          
          if (filter(idx1d)) {
            pos[e1_dim] = (0)*cell_size[e1_dim];
            pos[e2_dim] = (e2 + 0.5)*cell_size[e2_dim];
            pos[e3_dim] = (e3 + 0.5)*cell_size[e3_dim];

            // points->InsertNextPoint(pos);
            //
            // orient_array->InsertNextTuple(orient_vec_opp);
            // post(idx1d);
          }
          
          for (; e1 < e1_count - 1; ++e1) {
            idx1d = ijk.dot(map_idx);
            bool idx1d_filter = filter(idx1d);

            auto adj_idx1d = idx1d + adj_step;
            bool adj_idx1d_filter = filter(adj_idx1d);

            if (idx1d_filter != adj_idx1d_filter) {
              pos[e1_dim] = (e1 + 1)*cell_size[e1_dim];
              pos[e2_dim] = (e2 + 0.5)*cell_size[e2_dim];
              pos[e3_dim] = (e3 + 0.5)*cell_size[e3_dim];

              

              if (idx1d_filter) {
                points->InsertNextPoint(pos);
                
                orient_array->InsertNextTuple(orient_vec);
                post(idx1d);
              }
              else {
                // points->InsertNextPoint(pos);
                //
                // // TODO
                // orient_array->InsertNextTuple(orient_vec_opp);
                // post(adj_idx1d);
              }
            }
          }


          idx1d = ijk.dot(map_idx);
          if (filter(idx1d)) {
            pos[e1_dim] = (e1_count)*cell_size[e1_dim];
            pos[e2_dim] = (e2 + 0.5)*cell_size[e2_dim];
            pos[e3_dim] = (e3 + 0.5)*cell_size[e3_dim];

            points->InsertNextPoint(pos);

            orient_array->InsertNextTuple(orient_vec);
            post(idx1d);
          }
        }
    }

    



    
  public:
    void LoadImage() {
      auto image_path = 
        // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
        // R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)"
        R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\EstNorm_500x500x500_4p0um.raw)"
      ;

      auto velems_path = 
        // R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
        R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
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
      vtkNew<vtkIntArray> velem_array;

      vtkNew<vtkIntArray> velem_adjacent_array;

      
      
      size_t size;
      dpl::vector_n<pnm_idx, 3> dim;
      
      
      {
        // Open the stream
        std::ifstream is(image_path);
        // Determine the file length
        is.seekg(0, std::ios_base::end);
        size = is.tellg();
        is.seekg(0, std::ios_base::beg);
        std::vector<unsigned char> v(size);
        is.read(reinterpret_cast<char*>(v.data()), size);


        
        for (size_t i = 0; i < size; ++i)
          phase_array->InsertTypedComponent(i, 0, v[i]/*v[i] == 0 ? 0 : 1*/);

        
          // phase_array->InsertTypedComponent(i, 0, v[i] == 0 ? 230 : 99);
          // phase_array->InsertNextTuple1(/*i<(size/2) ? 1 : 0*/v[i] == 0 ? 230 : 99);

        phase_array->SetName("phase");
        // image_data_->GetCellData()->SetScalars(phase_array);
        image_data_->GetCellData()->AddArray(phase_array);
      }

      dim = std::round(std::cbrt(size));

      {
        pore_network_model pnm;
        auto velems = pnm.read_icl_velems(velems_path, dim);


        std::map<int, int> group;
        for (pnm_idx i = 0, count = velems.size(); i < count; ++i) {
          auto val = velems[i];

          ++group[val];

          velem_array->InsertTypedComponent(i, 0, val);
        }    
        
        
        velem_array->SetName("velem");
        image_data_->GetCellData()->AddArray(velem_array);
        

        auto [min, max] = std::ranges::minmax_element(velems);
        auto count = *max + 1;

        
        lut_velem_->IndexedLookupOn();
        lut_velem_->SetNumberOfTableValues(count);
        for (int32_t i = 0; i < count; ++i) {
          auto coef = static_cast<double>(i*45%count)/count;

          auto color = QColor::fromHsl(coef*255, 175, 122);
          lut_velem_->SetTableValue(i, color.redF(), color.greenF(), color.blueF());  // NOLINT(clang-diagnostic-double-promotion)
          lut_velem_->SetAnnotation(vtkVariant(i), std::to_string(i));
        }

        lut_velem_->SetTableValue(0, 0.4, 0.4, 0.4);
        lut_velem_->SetTableValue(1, 0.55, 0.55, 0.55);
        

        
        // dpl::vtk::PopulateLutRedWhiteBlue(lut_velem_);
        // lut_velem_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);


        {
          pnm_3idx map_idx{1, dim.x(), dim.x()*dim.y()};
          pnm_3idx ijk;

          pnm_idx idx1d = 0;

          auto& i = ijk.x();
          auto& j = ijk.y();
          auto& k = ijk.z();

          for (k = 0; k < dim.z(); ++k)
            for (j = 0; j < dim.y(); ++j)
              for (i = 0; i < dim.x(); ++i, ++idx1d) {
                int32_t adj = 1;
                
                if (velems[idx1d] < 2) {
                  if (i > 0) {
                    auto adj_idx = (ijk - pnm_3idx{1, 0, 0}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }

                  if (i < dim.x() - 1) {
                    auto adj_idx = (ijk + pnm_3idx{1, 0, 0}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }


                  if (j > 0) {
                    auto adj_idx = (ijk - pnm_3idx{0, 1, 0}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }

                  if (j < dim.y() - 1) {
                    auto adj_idx = (ijk + pnm_3idx{0, 1, 0}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }


                  if (k > 0) {
                    auto adj_idx = (ijk - pnm_3idx{0, 0, 1}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }

                  if (k < dim.z() - 1) {
                    auto adj_idx = (ijk + pnm_3idx{0, 0, 1}).dot(map_idx);
                    if (velems[adj_idx] > 1)
                      adj = velems[adj_idx];
                  }
                }
                else {
                  adj = 0;
                }

                velem_adjacent_array->InsertTypedComponent(idx1d, 0, adj);
                // ++idx;

                // auto val = file_ptr[velems_factor.dot(ijk + 1)];
                // *velems_ptr++ = val < 0 ? val + 2 : val;
              }


          velem_adjacent_array->SetName("velem_adj");
          image_data_->GetCellData()->AddArray(velem_adjacent_array);
        }
            
      }



     

      

      



      

      
      
      
      auto dim_side = std::round(std::cbrt(size));
      image_data_->SetDimensions(v3i{dim_side + 1});
      image_data_->SetSpacing(v3d{1.0/dim_side});
      image_data_->SetOrigin(v3d{0});

      
      


      vtkDataArray* selected_arr = velem_adjacent_array;

      // vtkDataArray* selected_arr = velem_array;
      
      
      image_data_->GetCellData()->SetActiveScalars(selected_arr->GetName());

      vtkNew<vtkDataSetMapper> mapper;
      
      {
        mapper->SetInputData(image_data_);
        mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, selected_arr->GetName()); // This is needed only for vtkImageData without vtkThreshold
      
        mapper->SetColorModeToMapScalars(); 
        mapper->UseLookupTableScalarRangeOn();
        mapper->SetScalarModeToUseCellData();
        mapper->SetLookupTable(image_data_->GetCellData()->GetScalars() == phase_array ? lut_pore_solid_ : lut_velem_);
      }
      
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
      



      
      
      

      
      
      image_actor_->SetMapper(mapper);
      image_actor_->GetProperty()->SetEdgeVisibility(false/*true*/);
      image_actor_->GetProperty()->SetEdgeColor(v3d{0.5} /*0, 0, 0*/);
      
      image_actor_->GetProperty()->SetAmbient(0.5);
      image_actor_->GetProperty()->SetDiffuse(0.4);
      image_actor_->GetProperty()->BackfaceCullingOn();
      

      // renderer_->AddActor(image_actor_);
    }

    
    void Init() {
      auto pnm_path = 
        // R"(E:\hwu\research126\d\modelling\networks\3D_network\10x10x10\10x10x10)"
        // R"(C:\dev\.temp\images\SS-1000\XNet)"
        // R"(E:\hwu\research126\d\modelling\networks\TwoScaleNet\MulNet)"
        // R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
        R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)"
      ;

      // v3i dim = 252;
      v3i dim = 500;

      
      qvtk_widget_ = new QVTKWidgetRef;
      qvtk_widget_->setRenderWindow(render_window_);
      render_window_->AddRenderer(renderer_);
      interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(v3d{1});


      // {
      //   tree_view_ = new QPropertyTreeView;
      //
      //   auto* hsplit = new QSplitter{Qt::Horizontal};
      //   hsplit->addWidget(tree_view_);
      //   hsplit->addWidget(qvtk_widget_);
      //   hsplit->setStretchFactor(0, 0);
      //   hsplit->setStretchFactor(1, 1);
      //
      //   setCentralWidget(hsplit);
      // }


      setCentralWidget(qvtk_widget_);


      

     



      


      

      



      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      
      // // pni_.read_from_binary_file(
      // //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\_MY_TEST_FILE_2.bin)"
      // //   // R"(C:\dev\.temp\_MY_TEST_FILE_2.bin)"
      // // );
      // // pni_.physical_size = {252};

      pore_network_model icl_pnm_inv{pnm_path, pore_network_model::file_format::statoil};























      
      
      auto [min, max] = std::ranges::minmax_element(icl_pnm_inv.node_.range(attribs::r_ins));
      
      lut_continuous_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);


      
      renderer_->AddActor(CreateNetworkAssembly(icl_pnm_inv, lut_continuous_));

      
      
      
      LoadImage();

      
      // image_actor_->SetUserTransform()
      {
        vtkNew<vtkTransform> trans;
        trans->PostMultiply();
        trans->Scale(v3d{icl_pnm_inv.physical_size.x()});
        // trans->Translate(icl_pnm_inv.physical_size.x(), 0, 0);
        image_actor_->SetUserTransform(trans);
      }



      
     

     

      
      

    

      
                // {
                //
                //   vtkNew<vtkPoints> points;
                //   vtkNew<vtkCellArray> cells;
                //   vtkNew<vtkPolyData> poly;
                //   poly->SetPoints(points);
                //   poly->SetPolys(cells);
                //
                //
                //   
                //   
                //
                //   auto* phase_in = static_cast<vtkUnsignedCharArray*>(image_data_->GetCellData()->GetArray("phase"));
                //   auto* velems_arr_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("velem"));
                //   auto* velems_adj_arr_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("velem_adj"));
                //
                //   vtkNew<vtkIntArray> velems_arr_out;
                //   velems_arr_out->SetName("velem_adj");
                //   
                //
                //   poly->GetCellData()->AddArray(velems_arr_out);
                //   
                //   
                //   auto pred = [&, this](pnm_idx idx) {
                //
                //
                //     // return true;
                //     
                //     // return phase_in->GetTypedComponent(idx, 0) == 2;
                //
                //
                //     
                //     return velems_arr_in->GetTypedComponent(idx, 0) < 2;
                //     // ;
                //     //
                //
                //     int z = idx/(dim.x()*dim.y());
                //     int y = (idx - z*dim.x()*dim.y())/dim.x();
                //     int x = idx - z*dim.x()*dim.y() - y*dim.x();
                //     
                //     // int z = i%dim.z();
                //     // int y = (i/dim.z())%dim.y();
                //     // int x = i/(dim.y()*dim.z()); 
                //
                //     
                //     
                //     return x < 250 && idx%2;
                //     
                //
                //
                //     // return i < ((252*252*252)/2);
                //   };
                //
                //   auto post = [&, this](pnm_idx idx) {
                //     auto val = velems_adj_arr_in->GetTypedComponent(idx, 0);
                //     velems_arr_out->InsertNextTypedTuple(&val);
                //   };
                //   
                //   FilterFaces<0>(dim, icl_pnm_inv.physical_size/dim, pred, post, poly);
                //   FilterFaces<1>(dim, icl_pnm_inv.physical_size/dim, pred, post, poly);
                //   FilterFaces<2>(dim, icl_pnm_inv.physical_size/dim, pred, post, poly);
                //
                //
                //   vtkNew<vtkPolyDataMapper> mapper;
                //   // mapper->SetInputData(poly);
                //   // mapper->SetScalarRange(poly->GetScalarRange());
                //
                //
                //   mapper->SetInputData(poly);
                //   // mapper->SetScalarRange(poly->GetScalarRange());
                //   mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, velems_arr_out->GetName()); // This is needed only for vtkImageData without vtkThreshold
                //   
                //   mapper->SetColorModeToMapScalars(); 
                //   mapper->UseLookupTableScalarRangeOn();
                //   mapper->SetScalarModeToUseCellData();
                //   mapper->SetLookupTable(/*image_data_->GetCellData()->GetScalars() == phase_array ? lut_pore_solid_ : */lut_velem_);
                //
                //   poly->GetCellData()->SetActiveScalars(velems_arr_out->GetName());
                //
                //
                //   vtkNew<vtkActor> actor;
                //   actor->SetMapper(mapper);
                //
                //   actor->GetProperty()->SetEdgeVisibility(false/*true*/);
                //   actor->GetProperty()->SetEdgeColor(v3d{0.5} /*0, 0, 0*/);
                //   
                //   actor->GetProperty()->SetAmbient(0.5);
                //   actor->GetProperty()->SetDiffuse(0.4);
                //   actor->GetProperty()->BackfaceCullingOn();
                //
                //   
                //   // renderer_->AddActor(actor);
                // }

      



      



      {

        
        vtkNew<vtkPolyData> quad;

        {
          vtkNew<vtkPoints> points;
          vtkNew<vtkCellArray> cells;
          quad->SetPoints(points);
          quad->SetPolys(cells);
          
          v3d pos;
          pos[0] = 0;

          pos[1] = -0.5;
          pos[2] = -0.5;
          points->InsertNextPoint(pos);

          pos[1] = 0.5;
          points->InsertNextPoint(pos);

          pos[2] = 0.5;
          points->InsertNextPoint(pos);

          pos[1] = -0.5;
          points->InsertNextPoint(pos);

          cells->InsertNextCell({0, 1, 2, 3});
        }

        


        {
          vtkNew<vtkPolyData> polydata;


          

          vtkNew<vtkGlyph3DMapper> glyphs;
          glyphs->SetSourceData(quad);
          glyphs->SetInputData(polydata);

         
          vtkNew<vtkFloatArray> orient_array;
          orient_array->SetName("orient");
          orient_array->SetNumberOfComponents(3);
          
          polydata->GetPointData()->AddArray(orient_array);
          glyphs->SetOrientationArray(orient_array->GetName());
          glyphs->SetOrientationModeToDirection();

          glyphs->SetScaleFactor(icl_pnm_inv.physical_size.x()/dim.x());
          

          vtkNew<vtkPoints> points;



          {
            auto* phase_in = static_cast<vtkUnsignedCharArray*>(image_data_->GetCellData()->GetArray("phase"));
            auto* velems_arr_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("velem"));
            auto* velems_adj_arr_in = static_cast<vtkIntArray*>(image_data_->GetCellData()->GetArray("velem_adj"));

            vtkNew<vtkIntArray> velems_arr_out;
            velems_arr_out->SetName("velem_adj");
        
            // poly->GetCellData()->AddArray(velems_arr_out);
        
        
            auto pred = [&, this](pnm_idx idx) {


              // return true;
              
              // return phase_in->GetTypedComponent(idx, 0) == 2;


              
              return velems_arr_in->GetTypedComponent(idx, 0) < 2;
              // ;/
              //

              int z = idx/(dim.x()*dim.y());
              int y = (idx - z*dim.x()*dim.y())/dim.x();
              int x = idx - z*dim.x()*dim.y() - y*dim.x();
              
              // int z = i%dim.z();
              // int y = (i/dim.z())%dim.y();
              // int x = i/(dim.y()*dim.z()); 

              
              
              return x < 250 && idx%2;
              


              // return i < ((252*252*252)/2);
            };

            auto post = [&, this](pnm_idx idx) {
              auto val = velems_adj_arr_in->GetTypedComponent(idx, 0);
              velems_arr_out->InsertNextTypedTuple(&val);
            };
            
            FilterFacesGlyph<0>(dim, icl_pnm_inv.physical_size/dim, pred, post, orient_array, points);
            // FilterFacesGlyph<1>(dim, icl_pnm_inv.physical_size/dim, pred, post, orient_array, points);
            // FilterFacesGlyph<2>(dim, icl_pnm_inv.physical_size/dim, pred, post, orient_array, points);





             polydata->GetPointData()->SetScalars(velems_arr_out);




            glyphs->SetLookupTable(lut_velem_);
            glyphs->SetColorModeToMapScalars();
            glyphs->UseLookupTableScalarRangeOn();
            glyphs->SetScalarModeToUsePointData();

            
            
          }





          // vtkNew<vtkFloatArray> color_array;
          // color_array->SetName("color");
          // color_array->SetNumberOfComponents(1);
         
          

          

          



          
          // orient_array->InsertNextTuple(v3d{0, 0, 1});
          // orient_array->InsertNextTuple(v3d{1, 0, 1});

          
          // points->InsertNextPoint(v3d{0, 0, 0});
          // points->InsertNextPoint(v3d{0, 0, 0});
          

          polydata->SetPoints(points);

          
          vtkNew<vtkActor> actor;
          actor->SetMapper(glyphs);
            

          actor->PickableOff();
          actor->DragableOff();
          actor->DebugOff();
          actor->GlobalWarningDisplayOff();
          

          actor->GetProperty()->SetEdgeVisibility(/*false*/true);
          actor->GetProperty()->SetEdgeColor(v3d{0.25} /*0, 0, 0*/);
          
          actor->GetProperty()->SetAmbient(0.5);
          actor->GetProperty()->SetDiffuse(0.4);
          actor->GetProperty()->BackfaceCullingOff();


          
          // actor->buffer();
          
          
          // renderer_->AddActor(actor);

          std::cout << "PRE UPD\n";

          
          
          glyphs->Update();
          // glyphs->MaskingOff();
          glyphs->OrientOff();
          
          actor->DragableOff();
          actor->DebugOff();
          actor->PickableOff();
          actor->UseBoundsOff();

          render_window_->Initialize();
          
          // renderer_->Render();
          
          std::cout << "POST UPD\n";
          
        }
        
      }












      








      

      
      
      

      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());
      tidy_axes_.SetFormat(".2e");

      

      // double bounds[] = {
      //   0., icl_pnm_inv.physical_size.x(),
      //   0., icl_pnm_inv.physical_size.y(),
      //   0., icl_pnm_inv.physical_size.z()};
      // tidy_axes_.Build(bounds);

      tidy_axes_.Build();
      
      // renderer_->ResetCamera();
    }
  };
}
