#pragma once

#include "xpm/functions.h"

#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>

// #include <QWidget>
#include <QMainWindow>
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

    vtkNew<vtkRenderer> renderer_; 
    vtkNew<vtkGenericOpenGLRenderWindow> render_window_;
    
    vtkRenderWindowInteractor* interactor_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;
    vtkNew<vtkLookupTable> lut_continuous_;
    vtkNew<vtkLookupTable> lut_pore_solid_;
    vtkNew<vtkLookupTable> lut_velem_;

    
    
    
    vtkNew<vtkImageData> image_data_;
    vtkNew<vtkThreshold> threshold_;



    vtkNew<vtkActor> image_actor_;




    


    
    

    


    static auto CreateNetworkAssembly(const pore_network_model& pnm, vtkLookupTable* lut) {
      auto net = vtkSmartPointer<vtkAssembly>::New();
      net->AddPart(CreateNodeActor(pnm, lut));
      net->AddPart(CreateThroatActor(pnm, lut));
      return net;
    }
    
  public:
    void LoadImage() {
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
        std::ifstream is(
          // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
          R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\images\Bmps252_6um.raw)"
        );
        // Determine the file length
        is.seekg(0, std::ios_base::end);
        size = is.tellg();
        is.seekg(0, std::ios_base::beg);
        std::vector<unsigned char> v(size);
        is.read(reinterpret_cast<char*>(v.data()), size);


        
        for (size_t i = 0; i < size; ++i)
          phase_array->InsertTypedComponent(i, 0, v[i] == 0 ? 0 : 1);
          // phase_array->InsertTypedComponent(i, 0, v[i] == 0 ? 230 : 99);
          // phase_array->InsertNextTuple1(/*i<(size/2) ? 1 : 0*/v[i] == 0 ? 230 : 99);

        phase_array->SetName("phase");
        // image_data_->GetCellData()->SetScalars(phase_array);
        image_data_->GetCellData()->AddArray(phase_array);
      }

      dim = std::round(std::cbrt(size));

      {
        pore_network_model pnm;
        auto velems = pnm.read_icl_velems(R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)", dim);


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
        lut_velem_->SetTableValue(1, 0.7, 0.7, 0.7);
        

        
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



      threshold_->SetInputData(image_data_);
      threshold_->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, velem_adjacent_array->GetName());
      
      // threshold_->SetLowerThreshold(2);
      // threshold_->SetUpperThreshold(100000);
      
      threshold_->SetLowerThreshold(1);
      threshold_->SetUpperThreshold(1e9);



      vtkNew<vtkDataSetMapper> mapper;

      mapper->SetInputData(threshold_->GetOutput());
      
      // mapper->SetInputData(image_data_);
      // mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, selected_arr->GetName()); // This is needed only for vtkImageData without vtkThreshold

      mapper->SetColorModeToMapScalars(); 
      mapper->UseLookupTableScalarRangeOn();
      mapper->SetScalarModeToUseCellData();
      mapper->SetLookupTable(image_data_->GetCellData()->GetScalars() == phase_array ? lut_pore_solid_ : lut_velem_);



      image_data_->Modified();
      threshold_->Modified();
      threshold_->Update();
      
      

      
      
      image_actor_->SetMapper(mapper);
      image_actor_->GetProperty()->SetEdgeVisibility(true);
      image_actor_->GetProperty()->SetEdgeColor(v3d{0.5} /*0, 0, 0*/);
      
      image_actor_->GetProperty()->SetAmbient(0.5);
      image_actor_->GetProperty()->SetDiffuse(0.4);
      image_actor_->GetProperty()->BackfaceCullingOn();
      

      renderer_->AddActor(image_actor_);
    }

    
    void Init() {
      qvtk_widget_ = new QVTKWidgetRef;
      qvtk_widget_->setRenderWindow(render_window_);
      render_window_->AddRenderer(renderer_);
      interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(v3d{1});

      this->setCentralWidget(qvtk_widget_);



      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      
      // // pni_.read_from_binary_file(
      // //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\_MY_TEST_FILE_2.bin)"
      // //   // R"(C:\dev\.temp\_MY_TEST_FILE_2.bin)"
      // // );
      // // pni_.physical_size = {252};

      pore_network_model icl_pnm_inv{
        // R"(E:\hwu\research126\d\modelling\networks\3D_network\10x10x10\10x10x10)"
        // R"(C:\dev\.temp\images\SS-1000\XNet)"
        // R"(E:\hwu\research126\d\modelling\networks\TwoScaleNet\MulNet)"
        R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)"
      
        , pore_network_model::file_format::statoil
      };

      
      
      auto [min, max] = std::ranges::minmax_element(icl_pnm_inv.node_.range(attribs::r_ins));
      
      lut_continuous_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);


      
      renderer_->AddActor(CreateNetworkAssembly(icl_pnm_inv, lut_continuous_));

      
      // pore_network_model icl_pnm{
      //   R"(C:\dev\pnextract\out\build\x64-Release\IMAGE1\)",
      //   pore_network_model::file_format::statoil
      // };
      
      // {
      //   auto vtk_pnm = CreateNetworkAssembly(icl_pnm, lut_continuous_);
      //   vtkNew<vtkTransform> trans;
      //   trans->PostMultiply();
      //   // trans->Scale(v3d{icl_pnm_inv.physical_size.x()});
      //   trans->Translate(2*icl_pnm_inv.physical_size.x(), 0, 0);
      //   vtk_pnm->SetUserTransform(trans);
      //   renderer_->AddActor(vtk_pnm);
      // }



                
                // pore_network_model uow_pnm{
                //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\_MY_TEST_FILE.bin)",
                //   pore_network_model::file_format::binary_dp69
                // };
                // uow_pnm.scale(icl_pnm_inv.physical_size.x()/252.0);
                //
                // {
                //   auto vtk_pnm = CreateNetworkAssembly(uow_pnm, lut_continuous_);
                //   vtkNew<vtkTransform> trans;
                //   trans->PostMultiply();
                //   trans->Translate(0, -1.1*icl_pnm_inv.physical_size.y(), 0);
                //   vtk_pnm->SetUserTransform(trans);
                //   renderer_->AddActor(vtk_pnm);
                // }

      
                //
                // pore_network_model uow_pnm_inv{
                //   R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\_MY_TEST_FILE_INV.bin)",
                //   pore_network_model::file_format::binary_dp69
                // };
                // uow_pnm_inv.scale(icl_pnm_inv.physical_size.x()/252.0);
                //
                // {
                //   auto vtk_pnm = CreateNetworkAssembly(uow_pnm_inv, lut_continuous_);
                //   vtkNew<vtkTransform> trans;
                //   trans->PostMultiply();
                //   trans->Translate(2*icl_pnm_inv.physical_size.x(), -1.1*icl_pnm_inv.physical_size.y(), 0);
                //   vtk_pnm->SetUserTransform(trans);
                //   renderer_->AddActor(vtk_pnm);
                // }
      
      
      LoadImage();

      
      // image_actor_->SetUserTransform()
      {
        vtkNew<vtkTransform> trans;
        trans->PostMultiply();
        trans->Scale(v3d{icl_pnm_inv.physical_size.x()});
        // trans->Translate(icl_pnm_inv.physical_size.x(), 0, 0);
        image_actor_->SetUserTransform(trans);
      }
      
      

      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());
      tidy_axes_.SetFormat(".2e");
      tidy_axes_.Build();

      // renderer_->ResetCamera();
    }
  };
}
