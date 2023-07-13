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
    
    // char processor_name[MPI_MAX_PROCESSOR_NAME];
    // int name_len;
    // MPI_Get_processor_name(processor_name, &name_len);

    // std::cout << std::format("Hello world from {}/{} (proc ID {})\n",
    //   mpi_rank, mpi_size, GetCurrentProcessId()) << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);


    // using namespace std::chrono;
    // using tp = time_point<steady_clock>;
    // tp t0, t1, t2, t3;

    // if (w_rank == root)
    //   t0 = std::chrono::high_resolution_clock::now();
    
    // dpl::hypre::InputDeprec input;
    dpl::hypre::ls_known_ref lk_ref;

    shared_memory_object smo{open_only, "xpm-hypre-input", read_only};
    mapped_region region(smo, read_only);


    std::pair<HYPRE_BigInt, HYPRE_BigInt>* range_ptr;

    dpl::hypre::load(region, lk_ref, range_ptr);
    


    // input.load(smo);

    // if (w_rank == root)
    //   t1 = std::chrono::high_resolution_clock::now();


    #ifdef HYPRE_SEQUENTIAL
      static constexpr auto jlower = 0;
      const auto jupper = lk_ref.nrows - 1;
    #else
      dpl::hypre::mpi_block::range = *range_ptr;
      auto [jlower, jupper] = dpl::hypre::mpi_block::range; //dpl::hypre::mpi_part(lk_ref.nrows);
    #endif


    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == root) {
      // std::cout << std::format("rank {}, range {}--{}", mpi_rank, dpl::hypre::mpi_block::range.first, dpl::hypre::mpi_block::range.second) << std::flush;
      // std::cout << "\n\nPRE_GATHER" << std::flush;
    // }
    // MPI_Barrier(MPI_COMM_WORLD);


    auto count = jupper - jlower + 1;
    
    auto indices = std::make_unique<HYPRE_BigInt[]>(count);
    for (HYPRE_BigInt i = 0; i < count; ++i)
      indices[i] = jlower + i;
    auto pressure_part = std::make_unique<HYPRE_Complex[]>(count);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == root)
    //   std::cout << "\n\nPRE_SOLVE" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    // try {

    dpl::hypre::solve(
      lk_ref/*.get_ref()*/,
      dpl::hypre::ls_unknown_ref{
        count,
        indices.get(),
        pressure_part.get()
      });

    // }
    // catch (std::exception& e){
    //   std::cout << std::format("exception in rank {}", mpi_rank) << std::flush;
    // }

    // MPI_Barrier(MPI_COMM_WORLD);


    // MPI_Barrier(MPI_COMM_WORLD);
    // if (w_rank == root) {
    //   std::cout << "\n\nPOST_SOLVE" << std::flush;
    // }
    // MPI_Barrier(MPI_COMM_WORLD);


    // auto pressure_part = input.Solve();

    // if (w_rank == root)
    //   t2 = std::chrono::high_resolution_clock::now();

    std::unique_ptr<double[]> pressure_utpr;
    // std::vector<double> pressure;
    double* recvbuf = nullptr; 


    std::unique_ptr<int[]> recvcounts_utpr;
    std::unique_ptr<int[]> displs_utpr;
    int* recvcounts = nullptr;
    int* displs = nullptr;
    
    if (mpi_rank == root) {
      recvbuf = new double[lk_ref.nrows];
      recvcounts = new int[mpi_size];
      displs = new int[mpi_size];

      pressure_utpr.reset(recvbuf);
      recvcounts_utpr.reset(recvcounts);
      displs_utpr.reset(displs);


      for (int i = 0; i < mpi_size; ++i) {
        recvcounts[i] = range_ptr[i].second - range_ptr[i].first + 1;
        displs[i] = range_ptr[i].first;
      }
    }

    
    
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == root)
    //   std::cout << "\n\nPRE_GATHER" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    
    

    MPI_Gatherv(
      pressure_part.get(), count, MPI_DOUBLE,
      recvbuf, recvcounts, displs, MPI_DOUBLE,
      root, MPI_COMM_WORLD);


    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == root)
    //   std::cout << "\n\nPOST_GATHER" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    // if (w_rank == root)
    //   t3 = std::chrono::high_resolution_clock::now();
  
    {
      // auto from_pnm = pnm.GenerateInput();
      
      if (mpi_rank == root) {
        // auto sum = std::accumulate(pressure.begin(), pressure.end(), 0.0);


        shared_memory_object smo{open_or_create, "xpm-hypre-output", read_write};
        smo.truncate(lk_ref.nrows*sizeof(double));
        mapped_region region(smo, read_write);
        std::memcpy(region.get_address(), recvbuf, lk_ref.nrows*sizeof(double));


        // std::cout <<
        //   std::format("\n\nLoad {}ms, Solve {}ms, MPI_Gather {}ms",
        //     duration_cast<milliseconds>(t1 - t0).count(),
        //     duration_cast<milliseconds>(t2 - t1).count(),
        //     duration_cast<milliseconds>(t3 - t2).count()
        //   );
        
        // shared_memory_object smo{open_or_create, "xpm-hypre-output", read_write};
        // smo.truncate(input.nrows*sizeof(double));
        // mapped_region region(smo, read_write);
        // auto* ptr = region.get_address();
        // std::memcpy()
        // *(double*)ptr = sum;
        
        // std::cout << sum << '\n';
      }
    }
    
    MPI_Finalize();

    // if (w_rank == root) {
    //   getchar();
    // }
    return 0;
  }
  else {
    MPI_Init(&argc, &argv);

    // {
    //
    //   xpm::pore_network_model pnm{
    //     R"(C:\dev\pnextract\out\build\x64-Release\Bmps252_INV\)",
    //     // R"(C:\dev\pnextract\out\build\x64-Release\EstThreePhase500_NORM\)",
    //     xpm::pore_network_model::file_format::statoil};
    //   
    //   auto input = pnm.GeneratePressureInput();
    //
    //
    //   // dpl::hypre::solve()
    //
    //   auto indices = std::make_unique<HYPRE_BigInt[]>(input.nrows);
    //   for (HYPRE_BigInt i = 0; i < input.nrows; ++i)
    //     indices[i] = i;
    //   std::vector<HYPRE_Real> vals(input.nrows);
    //
    //   dpl::hypre::solve(
    //     input.get_ref(),
    //     dpl::hypre::ls_unknown_ref{
    //       input.nrows,
    //       indices.get(),
    //       vals.data()
    //     });
    //
    //   // auto vals = input.Solve();
    //
    //   std::cout << std::format("\n\n {}", std::accumulate(vals.begin(), vals.end(), 0.0));
    //
    //   getchar();
    //
    // }


    // {
    //   shared_memory_object smo{open_or_create, "xpm-hypre-input", read_write};
    //   input.Save(smo);
    // }
    //
    //
    // auto solve_result = std::system(
    //   std::format("mpiexec -n 4 \"{}\" -s",
    //     "xpm_project.exe"
    //   ).c_str());
    //
    //
    // shared_memory_object smo_output{open_only, "xpm-hypre-output", read_only};
    // mapped_region mr_ouput{smo_output, read_only};
    // auto* ptr = static_cast<double*>(mr_ouput.get_address());
    //
    // auto val = std::accumulate(ptr, ptr + input.nrows, 0.0);


    
    

    // auto value = *(double*)mapped_region{
    //   shared_memory_object{open_only, "xpm-hypre-output", read_only}, read_only
    // }.get_address();

    // auto* ptr = region.get_address();

    
    


    
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