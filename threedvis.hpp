#pragma once

#include "dpl/vtk/TidyAxes.hpp"

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

    // vtkNew<vtkSphereSource> sphere_source_;
    // vtkNew<vtkActor> sphere_actor_;

    dpl::vtk::TidyAxes tidy_axes_;



    dpl::vector3d comp_angles(const dpl::vector3d& v) {
      auto norm = v.normalise();
      
      return {
        asin(norm.z())/3.141592654*180.0, 
        0.0f,                                                  
        (atan2(-norm.x(), norm.y())/3.141592654*180.0)      
      };
    }

    using v3d = dpl::vector3d;
    
  public:

    void Init() {
      qvtk_widget_ = new QVTKWidgetRef;
      qvtk_widget_->setRenderWindow(render_window_);
      render_window_->AddRenderer(renderer_);
      interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(v3d{1});

      this->setCentralWidget(qvtk_widget_);



      vtkNew<vtkPoints> points;
      points->InsertNextPoint(0, 0, 0.0);
      points->InsertNextPoint(3, -2, 0.0);
      points->InsertNextPoint(6, 2, -4);
      
      vtkNew<vtkPolyData> polydata;
      polydata->SetPoints(points);

      vtkNew<vtkConeSource> cylinder;
      cylinder->SetResolution(8);

      cylinder->SetCenter(0.5, 0.0, 0.0);
      cylinder->SetHeight(1);
      
      vtkNew<vtkGlyph3DMapper> glyph3DMapper;
      glyph3DMapper->SetSourceConnection(cylinder->GetOutputPort());
      glyph3DMapper->SetInputData(polydata);
      


      vtkNew<vtkFloatArray> orient_array;
      orient_array->SetName("orient");
      orient_array->SetNumberOfComponents(3);
      orient_array->InsertNextTuple(v3d{1, 0, 0});
      orient_array->InsertNextTuple(v3d{0, 1, 0});
      orient_array->InsertNextTuple(v3d{0, 0, -1});

      polydata->GetPointData()->AddArray(orient_array);
      glyph3DMapper->SetOrientationArray(orient_array->GetName());
      glyph3DMapper->SetOrientationModeToDirection();
      

      vtkNew<vtkFloatArray> scale_array;
      scale_array->SetName("scale");
      scale_array->SetNumberOfComponents(3);
      scale_array->InsertNextTuple(v3d{2, 1, 1});
      scale_array->InsertNextTuple(v3d{1, 2, 1});
      scale_array->InsertNextTuple(v3d{1, 2, 2});

      polydata->GetPointData()->AddArray(scale_array);
      glyph3DMapper->SetScaleArray(scale_array->GetName());
      glyph3DMapper->SetScaleModeToScaleByVectorComponents();

      
      
      
      

      vtkNew<vtkActor> actor;
      actor->SetMapper(glyph3DMapper);
      
      vtkNew<vtkNamedColors> colors;
      actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
      renderer_->AddActor(actor);
      

                  

      renderer_->ResetCamera();
      
      tidy_axes_.Init(renderer_.Get());
      tidy_axes_.Build();



      // // sphere_source_->SetRadius(3.0);
                  // // sphere_source_->SetPhiResolution(5);
                  // // sphere_source_->SetThetaResolution(5);
                  // // sphere_source_->Update();
                  //
                  //
                  // vtkNew<vtkCylinderSource> cylinder;
                  // cylinder->SetResolution(8);
                  //
                  // cylinder->SetCenter(0, 0.5, 0);
                  // cylinder->Update();
                  //
                  //
                  //
                  //
                  //
                  //
                  //
                  // // sphere_actor_->VisibilityOff();
                  //
                  //
                  //
                  // vtkNew<vtkPoints> points;
                  // points->InsertNextPoint(0, 0, 0.0);
                  // points->InsertNextPoint(7, 0, 0.0);
                  // points->InsertNextPoint(14, 0, 0.0);
                  //
                  // vtkNew<vtkPolyData> polydata;
                  // polydata->SetPoints(points);
                  //
                  //
                  //
                  //
                  // vtkNew<vtkGlyph3D> glyph3D;
                  // glyph3D->SetSourceConnection(cylinder->GetOutputPort());
                  // glyph3D->SetInputData(polydata);
                  //
                  // {
                  //   // vtkNew<vtkDoubleArray> rad_array;
                  //   // rad_array->InsertNextTuple1(1);
                  //   // rad_array->InsertNextTuple1(3);
                  //   // rad_array->InsertNextTuple1(5);
                  //   // polydata->GetPointData()->SetScalars(rad_array);
                  //   // glyph3D->SetScaleMode(VTK_SCALE_BY_SCALAR);
                  // }
                  //
                  //
                  // vtkNew<vtkFloatArray> rad_array;
                  // rad_array->SetNumberOfComponents(3);
                  // // rad_array->InsertNextTuple3(1, 10, 1);
                  // // rad_array->InsertNextTuple3(1, 5, 1);
                  // // rad_array->InsertNextTuple3(1, 20, 1);
                  //
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  // // polydata->GetPointData()->SetVectors(rad_array);
                  //
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  // // rad_array->InsertNextTuple3(0, 1, 0);
                  //
                  // rad_array->InsertNextTuple3(1, 1, 1);
                  // rad_array->InsertNextTuple3(1, 1, 1);
                  // rad_array->InsertNextTuple3(1, 1, 1);
                  // polydata->GetPointData()->SetNormals(rad_array);
                  // // polydata->GetPointData()->GetTCoords()
                  // // polydata->GetPointData()->normal
                  // // polydata->GetCellData()
                  //
                  // // glyph3D->SetScaleMode(VTK_SCALE_BY_VECTORCOMPONENTS);
                  // // glyph3D->SetScaleModeToScaleByVectorComponents();
                  // // glyph3D->SetVectorModeToUseVector();
                  // // glyph3D->SetVectorModeToUseNormal();
                  // // glyph3D->OrientOff();
                  // // glyph3D->ClampingOff();
                  //
                  // // glyph3D->SetVectorModeToUseVector();
                  //
                  //
                  // glyph3D->Update();
                  //
                  //
                  //
                  //
                  //
                  //
                  //
                  // vtkNew<vtkPolyDataMapper> mapper;
                  // mapper->SetInputConnection(glyph3D->GetOutputPort());
                  // vtkNew<vtkActor> actor;
                  // actor->SetMapper(mapper);
                  //
                  // vtkNew<vtkNamedColors> colors;
                  // actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
                  // renderer_->AddActor(actor);
      
      
      // vtkNew<vtkPolyDataMapper> mapper;
      // mapper->SetInputData(sphere_source_->GetOutput());
      // sphere_actor_->SetMapper(mapper);
      // sphere_actor_->GetProperty()->SetColor(0.5, 0.5, 0.5);
      // renderer_->AddActor(sphere_actor_);
    }
  };
}






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
