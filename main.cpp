// import threedvisMODULE;

#include "threedvis.hpp"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include <QWidget>
// #include <QStyleFactory>


#include <array>

// #include "ui_mywidg.h"

// struct my_t
// {
//   auto operator()() const {
//     return 22;
//   }
//
//   my_t() = default;
//   
//   my_t(const my_t& other) = delete;
//   my_t(my_t&& other) noexcept = delete;
//   my_t& operator=(const my_t& other) = delete;
//   my_t& operator=(my_t&& other) noexcept = delete;
// };
//
// auto complex(const auto& foo) {
//   return foo() + 3;
// }

int main(int argc, char* argv[])
{
  // auto qqq = complex([](){return 10;});
  //
  // auto my_val = my_t{};
  //
  // auto& my_val_ref = my_val;
  //
  // auto qqq2222 = complex(my_t{});
  //
  // auto qqq2222 = complex(my_val);
  //
  // auto qqq22 = complex(my_val_ref);
  
  
  
  // foo f;
  // f.helloworld();

  // constexpr auto qqq = dpl::vector3i{5, 6, 7};
  //
  // constexpr auto rmg = dpl::cdims<1>::tie(qqq);
  //
  // auto& [a, b, c] = rmg;
  //
  // std::array<double, std::get<1>(rmg)> qqewe;
  

  
  

  auto format = xpm::QVTKWidgetRef::defaultFormat();

  #ifdef _WIN32
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
  #else
    format.setProfile(QSurfaceFormat::CoreProfile);
  #endif
  

  #if (VTK_MAJOR_VERSION == 8)
    QSurfaceFormat::setDefaultFormat(format);
  #elif (VTK_MAJOR_VERSION == 9)
  #endif
  
  
  QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling, false);
  
  QApplication app(argc, argv);

  
  
  // QWidget widget;
  xpm::XPMWidget widget;
  widget.Init();
  
  // Ui::MainWindow ui;
  // ui.setupUi(&widget);

  widget.resize(1400, 1000);
  
  widget.show();
  return app.exec();

  
  // vtkNew<vtkNamedColors> colors;
  //
  // // Set the background color.
  // std::array<unsigned char, 4> bkg{{26, 51, 102, 255}};
  // colors->SetColor("BkgColor", bkg.data());
  //
  // // This creates a polygonal cylinder model with eight circumferential facets
  // // (i.e, in practice an octagonal prism).
  // vtkNew<vtkCylinderSource> cylinder;
  // cylinder->SetResolution(8);
  //
  // // The mapper is responsible for pushing the geometry into the graphics
  // // library. It may also do color mapping, if scalars or other attributes are
  // // defined.
  // vtkNew<vtkPolyDataMapper> cylinderMapper;
  // cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
  //
  // // The actor is a grouping mechanism: besides the geometry (mapper), it
  // // also has a property, transformation matrix, and/or texture map.
  // // Here we set its color and rotate it around the X and Y axes.
  // vtkNew<vtkActor> cylinderActor;
  // cylinderActor->SetMapper(cylinderMapper);
  // cylinderActor->GetProperty()->SetColor(
  //     colors->GetColor4d("Tomato").GetData());
  // cylinderActor->RotateX(30.0);
  // cylinderActor->RotateY(-45.0);
  //
  // // The renderer generates the image
  // // which is then displayed on the render window.
  // // It can be thought of as a scene to which the actor is added
  // vtkNew<vtkRenderer> renderer;
  // renderer->AddActor(cylinderActor);
  // renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
  // // Zoom in a little by accessing the camera and invoking its "Zoom" method.
  // renderer->ResetCamera();
  // renderer->GetActiveCamera()->Zoom(1.5);
  //
  // // The render window is the actual GUI window
  // // that appears on the computer screen
  // vtkNew<vtkRenderWindow> renderWindow;
  // renderWindow->SetSize(300, 300);
  // renderWindow->AddRenderer(renderer);
  // renderWindow->SetWindowName("Cylinder");
  //
  // // The render window interactor captures mouse events
  // // and will perform appropriate camera or actor manipulation
  // // depending on the nature of the events.
  // vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  // renderWindowInteractor->SetRenderWindow(renderWindow);
  //
  // // This starts the event loop and as a side effect causes an initial render.
  // renderWindow->Render();
  // renderWindowInteractor->Start();
  //
  // return EXIT_SUCCESS;
}
