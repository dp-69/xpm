// import threedvisMODULE;

#include "xpm_widget.hpp"

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

// #include <mpi.h>



int main(int argc, char* argv[])
{
  // MPI_Init(&argc, &argv);
  //
  //
  //
  // // Get the number of processes
  // int world_size;
  // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  //
  // // Get the rank of the process
  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  //
  //
  //
  //
  //
  // // Get the name of the processor
  // char processor_name[MPI_MAX_PROCESSOR_NAME];
  // int name_len;
  // MPI_Get_processor_name(processor_name, &name_len);
  //
  // // Print off a hello world message
  // printf("Hello world from processor %s, rank %d out of %d processors, %d Proc ID\n",
  //        processor_name, world_rank, world_size, GetCurrentProcessId());
  //
  // std::cout << std::flush;
  //
  // if (world_rank == 0) {
  //   getchar();  
  // }
  //
  // MPI_Barrier(MPI_COMM_WORLD);
  //
  // // std::this_thread::sleep_for(std::chrono::seconds{3});
  //  
  //
  // xpm::pore_network_model pnm{
  //   R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)",
  //   // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)",
  //   xpm::pore_network_model::file_format::statoil};
  //
  //
  // {
  //   auto pressure = pnm.SolvePressure();
  //   
  //   auto sum = std::accumulate(pressure.begin(), pressure.end(), 0.0);
  //   std::cout << sum << '\n';
  // }
  //
  //
  // // Finalize the MPI environment.
  // MPI_Finalize();
  //
  //
  //
  //
  // return 0;


  
  

  
  

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
