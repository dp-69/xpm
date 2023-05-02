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
#include <vtkAssembly.h>

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
      // lut_discrete_->SetTableRange(0, 1);
      
      lut_pore_solid_->Modified();


      vtkNew<vtkUnsignedCharArray> phase_array;
      vtkNew<vtkIntArray> velem_array;
      
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
          phase_array->InsertNextTuple1(/*i<(size/2) ? 1 : 0*/v[i] == 0 ? 0 : 1);
        phase_array->SetName("phase");
        image_data_->GetCellData()->SetScalars(phase_array);
      }

      dim = std::round(std::cbrt(size));

      {
        pore_network_model pnm;
        auto velems = pnm.read_icl_velems(R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)", dim);


        std::map<int, int> group;
        for (pnm_idx i = 0, count = velems.size(); i < count; ++i) {
          ++group[velems[i]];
        }    
        

        int i = 0;
        for (auto& x : velems) {
          // x = ++i;
          velem_array->InsertNextTuple1(x);
        }
        velem_array->SetName("velem");
        image_data_->GetCellData()->AddArray(velem_array);
        

        auto [min, max] = std::ranges::minmax_element(velems);

        dpl::vtk::PopulateLutRedWhiteBlue(lut_velem_);
        lut_velem_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);

        
        
            
      }





      
      
      
      auto dim_side = std::round(std::cbrt(size));
      image_data_->SetDimensions(v3i{dim_side + 1});
      image_data_->SetSpacing(v3d{1.0/dim_side});
      image_data_->SetOrigin(v3d{0});

      
      

      image_data_->GetCellData()->SetActiveScalars(velem_array->GetName());


      

      vtkNew<vtkDataSetMapper> mapper;

      mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, velem_array->GetName());
      mapper->SetColorModeToMapScalars();
      // mapper->ColorByArrayComponent(vtkarr->GetName(), 0);
      // mapper->colo
      
      mapper->SetInputData(image_data_);
      mapper->SetScalarModeToUseCellData();
      // mapper->SetScalarModeToUseCellData();
      mapper->UseLookupTableScalarRangeOn();
      mapper->SetLookupTable(lut_velem_);

      // image_data_->Modified();
      // mapper->Modified();
      // mapper->Update();
      

      
      
      image_actor_->SetMapper(mapper);
      image_actor_->GetProperty()->SetEdgeVisibility(false);
      image_actor_->GetProperty()->SetEdgeColor(0, 0, 0);
      
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



                
                pore_network_model uow_pnm{
                  R"(C:\Users\dmytr\OneDrive - Heriot-Watt University\temp\_MY_TEST_FILE.bin)",
                  pore_network_model::file_format::binary_dp69
                };
                uow_pnm.scale(icl_pnm_inv.physical_size.x()/252.0);
                
                {
                  auto vtk_pnm = CreateNetworkAssembly(uow_pnm, lut_continuous_);
                  vtkNew<vtkTransform> trans;
                  trans->PostMultiply();
                  trans->Translate(0, -1.1*icl_pnm_inv.physical_size.y(), 0);
                  vtk_pnm->SetUserTransform(trans);
                  renderer_->AddActor(vtk_pnm);
                }
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
        trans->Translate(icl_pnm_inv.physical_size.x(), 0, 0);
        image_actor_->SetUserTransform(trans);
      }
      
      

      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());
      tidy_axes_.Build();

      // renderer_->ResetCamera();
    }
  };
}
