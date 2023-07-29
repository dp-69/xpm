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

#include <mpi.h>

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
// #include <boost/interprocess/windows_shared_memory.hpp>


int main(int argc, char* argv[])
{
  if (argc == 2 && !std::strcmp(argv[1], "-s")) {
    using namespace boost::interprocess;
    
    MPI_Init(&argc, &argv);
  
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    static constexpr auto root = 0;
    


    auto input = dpl::hypre::mpi::load(mpi_rank);


    // dpl::hypre::ls_known_ref lk_ref;
    // mapped_region region{shared_memory_object{open_only, dpl::hypre::mpi::smo_name, read_only}, read_only};
    // std::pair<HYPRE_BigInt, HYPRE_BigInt>* range_ptr;

    // auto nrows = dpl::hypre::load(region, lk_ref, range_ptr);

    // #ifdef HYPRE_SEQUENTIAL
    //   static constexpr auto jlower = 0;
    //   const auto jupper = lk_ref.nrows - 1;
    // #else
    //   dpl::hypre::mpi::range = *range_ptr;
    //   auto [jlower, jupper] = dpl::hypre::mpi::range; //dpl::hypre::mpi_part(lk_ref.nrows);
    // #endif

    auto [jlower, jupper] = *input.local_rows;

    auto local_nrows = jupper - jlower + 1;


    auto values = std::make_unique<double[]>(local_nrows);
    dpl::hypre::solve(jlower, jupper, input, values.get());


    std::unique_ptr<double[]> pressure_utpr;
    double* recvbuf = nullptr; 

    std::unique_ptr<int[]> recvcounts_utpr;
    std::unique_ptr<int[]> displs_utpr;
    int* recvcounts = nullptr;
    int* displs = nullptr;
    
    if (mpi_rank == root) {
      recvbuf = new double[input.global_nrows];
      recvcounts = new int[mpi_size];
      displs = new int[mpi_size];

      pressure_utpr.reset(recvbuf);
      recvcounts_utpr.reset(recvcounts);
      displs_utpr.reset(displs);


      for (int i = 0; i < mpi_size; ++i) {
        recvcounts[i] = input.local_rows[i].second - input.local_rows[i].first + 1;
        displs[i] = input.local_rows[i].first;
      }
    }

    MPI_Gatherv(
      values.get(), local_nrows, MPI_DOUBLE,
      recvbuf, recvcounts, displs, MPI_DOUBLE,
      root, MPI_COMM_WORLD);

    if (mpi_rank == root) {
      shared_memory_object smo_output{open_or_create, "xpm-hypre-output", read_write};
      smo_output.truncate(input.global_nrows*sizeof(double));
      mapped_region region_output(smo_output, read_write);
      std::memcpy(region_output.get_address(), recvbuf, input.global_nrows*sizeof(double));
    }

    MPI_Finalize();

    return 0;
  }
  else {
    MPI_Init(&argc, &argv);

    


    
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
    

    // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling, false);
    
    QApplication app(argc, argv);

    
    
    // QWidget widget;
    xpm::XPMWidget widget;

    
    
    widget.Init();
    
    // Ui::MainWindow ui;
    // ui.setupUi(&widget);

    widget.resize(1400, 1000);
    widget.show();

    auto result = QApplication::exec();

    MPI_Finalize();

    return result;
  }
}




// shared_memory_object shm (open_only, "MySharedMemory", read_only);
    // mapped_region region(shm, read_only);
    //
    // auto shared_inp = *((int*)region.get_address());
    //
    // // //Check that memory was initialized to 1
    // // char *mem = static_cast<char*>(region.get_address());
    // // for(std::size_t i = 0; i < region.get_size(); ++i)
    // //    if(*mem++ != 1)
    // //       return 1;   //Error checking memory
    //
    // return 5 + shared_inp;


  

  
  // struct shm_remove
  // {
  //    shm_remove() { shared_memory_object::remove("MySharedMemory"); }
  //    ~shm_remove(){ shared_memory_object::remove("MySharedMemory"); }
  // } remover;

  
  // shm.truncate(sizeof(int));
  // mapped_region region(shm, read_write);
  // *((int*)region.get_address()) = 20;
  
  


  
  
  
  
  
  
  
  
  // if (world_rank == 0) {
  //   getchar();  
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  
  // std::this_thread::sleep_for(std::chrono::seconds{3});




  
  
 


  


  
  

  
  



  

  

  
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