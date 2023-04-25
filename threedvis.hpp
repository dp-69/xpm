#pragma once

#include "xpm/pore_network_info.hpp"

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
    vtkNew<vtkLookupTable> lut_discrete_;
    vtkNew<vtkImageData> image_data_;

    vtkNew<vtkActor> image_actor_;

    pore_network_info pni_;


    static dpl::vector3d angles_for_j_norm(const dpl::vector3d& dir) {
      auto norm = dir.normalise();
      
      return {
        asin(norm.z())/3.141592654*180.0, 
        0.0f,                                                  
        (atan2(-norm.x(), norm.y())/3.141592654*180.0)      
      };
    }

    using v3i = dpl::vector3i;
    using v3d = dpl::vector3d;
    
  public:
    void LoadPoreNetwork() {
      pni_.read_from_binary_file(
        R"(C:\dev\.temp\_MY_TEST_FILE.bin)"
        // R"(C:\dev\.temp\_MY_TEST_FILE_2.bin)"
      );
      pni_.physical_size = {252};
      
      // pni_.read_from_text_file(
      //   // R"(E:\hwu\research126\d\modelling\networks\3D_network\10x10x10\10x10x10)"
      //   R"(C:\dev\.temp\images\SS-1000\XNet)"
      //   // R"(E:\hwu\research126\d\modelling\networks\TwoScaleNet\MulNet)"
      // );

      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);
      auto [min, max] = std::minmax_element(pni_.r_ins_.begin(), pni_.r_ins_.begin() + pni_.throat_count() + pni_.inner_node_count());
      lut_continuous_->SetTableRange(*min - (*max - *min)*0.1, *max + (*max - *min)*0.1);

      {
        vtkNew<vtkPolyData> polydata;
        
        vtkNew<vtkCylinderSource> cylinder;
        cylinder->SetResolution(20);
        cylinder->SetCenter(0.0, 0.5, 0.0);
        cylinder->SetRadius(1);
        cylinder->SetHeight(1);
        
        vtkNew<vtkGlyph3DMapper> throat_glyphs;
        throat_glyphs->SetSourceConnection(cylinder->GetOutputPort());
        throat_glyphs->SetInputData(polydata);
                       
        

        vtkNew<vtkFloatArray> orient_array;
        orient_array->SetName("orient");
        orient_array->SetNumberOfComponents(3);

        polydata->GetPointData()->AddArray(orient_array);
        throat_glyphs->SetOrientationArray(orient_array->GetName());
        throat_glyphs->SetOrientationModeToRotation();

        vtkNew<vtkFloatArray> scale_array;
        scale_array->SetName("scale");
        scale_array->SetNumberOfComponents(3);

        polydata->GetPointData()->AddArray(scale_array);
        throat_glyphs->SetScaleArray(scale_array->GetName());
        throat_glyphs->SetScaleModeToScaleByVectorComponents();



        
        

        vtkNew<vtkFloatArray> color_array;
        color_array->SetName("color");
        color_array->SetNumberOfComponents(1);
        polydata->GetPointData()->SetScalars(color_array);

        throat_glyphs->SetLookupTable(lut_continuous_);
        throat_glyphs->SetColorModeToMapScalars();
        throat_glyphs->UseLookupTableScalarRangeOn();
        throat_glyphs->SetScalarModeToUsePointData();

        vtkNew<vtkPoints> points;

        for (pnm_idx i = 0, count = pni_.throat_count(); i < count; ++i)
          if (auto [n0, n1] = pni_.throats[i];
            pni_.inner_node(n0) && pni_.inner_node(n1)) {
            auto& n0_pos = pni_.node_pos[n0];

            points->InsertNextPoint(n0_pos);
            orient_array->InsertNextTuple(angles_for_j_norm(pni_.node_pos[n1] - n0_pos));

            scale_array->InsertNextTuple(v3d{
              pni_.r_ins_[i], (pni_.node_pos[n1] - n0_pos).length(), pni_.r_ins_[i]
            });

            color_array->InsertNextTuple1(pni_.r_ins_[i]);
          }

        polydata->SetPoints(points);
        
        vtkNew<vtkActor> actor;
        actor->SetMapper(throat_glyphs);
        

        actor->GetProperty()->SetSpecularColor(1, 1, 1);
        actor->GetProperty()->SetSpecular(0);
        actor->GetProperty()->SetSpecularPower(75);
        actor->GetProperty()->SetAmbient(0.15);
        actor->GetProperty()->SetDiffuse(0.9);

        
        
        vtkNew<vtkNamedColors> colors;
        actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
        
        renderer_->AddActor(actor);
      }


      {
        vtkNew<vtkPolyData> polydata;
        
        vtkNew<vtkSphereSource> cylinder;
        cylinder->SetPhiResolution(10);
        cylinder->SetThetaResolution(10);
        cylinder->SetCenter(0, 0, 0);
        cylinder->SetRadius(1);
        
        vtkNew<vtkGlyph3DMapper> node_glyphs;
        node_glyphs->SetSourceConnection(cylinder->GetOutputPort());
        node_glyphs->SetInputData(polydata);
        node_glyphs->OrientOff();

        vtkNew<vtkFloatArray> scale_array;
        scale_array->SetName("scale");
        scale_array->SetNumberOfComponents(1);

        polydata->GetPointData()->AddArray(scale_array);
        node_glyphs->SetScaleArray(scale_array->GetName());
        node_glyphs->SetScaleModeToScaleByMagnitude();




        vtkNew<vtkFloatArray> color_array;
        color_array->SetName("color");
        color_array->SetNumberOfComponents(1);
        polydata->GetPointData()->SetScalars(color_array);

        node_glyphs->SetLookupTable(lut_continuous_);
        node_glyphs->SetColorModeToMapScalars();
        node_glyphs->UseLookupTableScalarRangeOn();
        node_glyphs->SetScalarModeToUsePointData();
        
        vtkNew<vtkPoints> points;
        
        for (pnm_idx i = 0, count = pni_.inner_node_count(); i < count; ++i) {
          points->InsertNextPoint(pni_.node_pos[i]);
          scale_array->InsertNextTuple1(pni_.node_r_ins(i));
          color_array->InsertNextTuple1(pni_.node_r_ins(i));
        }
        
        polydata->SetPoints(points);
        
        vtkNew<vtkActor> actor;
        actor->SetMapper(node_glyphs);

        actor->GetProperty()->SetSpecularColor(1, 1, 1);
        actor->GetProperty()->SetSpecular(0);
        actor->GetProperty()->SetSpecularPower(75);
        actor->GetProperty()->SetAmbient(0.15);
        actor->GetProperty()->SetDiffuse(0.9);
        
        vtkNew<vtkNamedColors> colors;
        actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
        renderer_->AddActor(actor);
      }
    }

    void LoadImage() {
      lut_discrete_->IndexedLookupOn();
      
      lut_discrete_->SetNumberOfTableValues(2);
      

      lut_discrete_->SetTableValue(0, 0.6, 0.6, 0.6);
      lut_discrete_->SetAnnotation(vtkVariant(0), "PORE");

      lut_discrete_->SetTableValue(1, 0, 0, 0);
      lut_discrete_->SetAnnotation(vtkVariant(1), "THROAT");
      // lut_discrete_->SetTableRange(0, 1);
      
      lut_discrete_->Modified();


      // Open the stream
      std::ifstream is(
        // R"(C:\dev\.temp\images\Bentheimer1000_normalized.raw)"
        R"(C:\dev\.temp\images\Bmps252_6um.raw)"
      );
      // Determine the file length
      is.seekg(0, std::ios_base::end);
      size_t size = is.tellg();
      is.seekg(0, std::ios_base::beg);
      std::vector<unsigned char> v(size);
      is.read(reinterpret_cast<char*>(v.data()), size);
      // Close the file
      is.close();


      vtkNew<vtkUnsignedCharArray> vtkarr;
      for (size_t i = 0; i < size; ++i)
        vtkarr->InsertNextTuple1(/*i<(size/2) ? 1 : 0*/v[i] == 0 ? 0 : 1);
      vtkarr->SetName("indicator");
      
      auto dim_side = std::round(std::cbrt(size));
      image_data_->SetDimensions(v3i{dim_side + 1});
      image_data_->SetSpacing(v3d{1.0/dim_side});
      image_data_->SetOrigin(v3d{0});

      image_data_->GetCellData()->SetScalars(vtkarr);
      

      

      vtkNew<vtkDataSetMapper> mapper;

      mapper->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkarr->GetName());
      mapper->SetColorModeToMapScalars();
      // mapper->ColorByArrayComponent(vtkarr->GetName(), 0);
      // mapper->colo
      
      mapper->SetInputData(image_data_);
      mapper->SetScalarModeToUseCellData();
      // mapper->SetScalarModeToUseCellData();
      mapper->UseLookupTableScalarRangeOn();
      mapper->SetLookupTable(lut_discrete_);

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




      
      LoadPoreNetwork();
      
      LoadImage();

      // image_actor_->SetUserTransform()
      vtkNew<vtkTransform> trans;
      trans->PostMultiply();
      trans->Scale(v3d{pni_.physical_size.x()});
      trans->Translate(pni_.physical_size.x(), 0, 0);
      image_actor_->SetUserTransform(trans);
      
      

      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());
      tidy_axes_.Build();

      // renderer_->ResetCamera();
    }
  };
}





// vtkNew<vtkPoints> points;
        // points->InsertNextPoint(0, 0, 0.0);
        // points->InsertNextPoint(3, -2, 0.0);
        // points->InsertNextPoint(6, 2, -4);
        
        

        // vtkNew<vtkCylinderSource> cylinder;
        // cylinder->SetResolution(8);
        //
        //
        // cylinder->SetCenter(0.0, 0.5, 0.0);
        // cylinder->SetHeight(1);
        //
        // vtkNew<vtkGlyph3DMapper> glyph3DMapper;
        // glyph3DMapper->SetSourceConnection(cylinder->GetOutputPort());
        // glyph3DMapper->SetInputData(polydata);
        //
        //
        //
        // vtkNew<vtkFloatArray> orient_array;
        // orient_array->SetName("orient");
        // orient_array->SetNumberOfComponents(3);
        //
        //
        //
        // orient_array->InsertNextTuple(angles_for_j_norm({1, 0, 0}));
        // orient_array->InsertNextTuple(angles_for_j_norm({0, 1, 0}));
        // orient_array->InsertNextTuple(angles_for_j_norm({0, 0, -1}));
        //
        // polydata->GetPointData()->AddArray(orient_array);
        // glyph3DMapper->SetOrientationArray(orient_array->GetName());
        // glyph3DMapper->SetOrientationModeToRotation();
        

        // vtkNew<vtkFloatArray> scale_array;
        // scale_array->SetName("scale");
        // scale_array->SetNumberOfComponents(3);
        // scale_array->InsertNextTuple(v3d{1, 1, 1});
        // scale_array->InsertNextTuple(v3d{1, 3, 1});
        // scale_array->InsertNextTuple(v3d{1, 2, 1});
        //
        // polydata->GetPointData()->AddArray(scale_array);
        // glyph3DMapper->SetScaleArray(scale_array->GetName());
        // glyph3DMapper->SetScaleModeToScaleByVectorComponents();






// // Global module fragment where #includes can happen
// module;
// #include <iostream>
//
// // first thing after the Global module fragment must be a module command
// export module fooMODULE;
//
// // char const* world() { return "hello world\n"; }
//
// export class foo {
// public:
//   foo();
//   ~foo();
//   void helloworld();
// };
//
// foo::foo() = default;
// foo::~foo() = default;
// void foo::helloworld() { std::cout << "hello world\n"; }
//
// // export module threedvisMODULE;
// //
// // #include <QWidget>;
// //
// //
// //
// //
// // namespace xpm
// // {
// //   export class FDWidget : public QWidget
// //   {
// //     
// //   };
// // }
