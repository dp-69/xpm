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

  std::vector<int> zomg;
  auto luuuq = zomg.end() - zomg.begin();

  
  xpm::voxel_idx pew{30};

  constexpr xpm::voxel_idx pew2{40};

  // constexpr bool qqqq = pew2 < pew;


  // auto qq = pew2 - pew;

  // xpm::voxel_idx omg = pew++;

  int k = 0;

  for (xpm::voxel_idx i{0}; *i < 10; ++i)
      k += *i;


  auto qqq = std::addressof(pew);
  auto qqqqq = std::addressof(++pew);


  auto powepw = static_cast<xpm::voxel_idx>(-1);


  // typename _Ty::_Signed_type;
                                           // typename _Ty::_Unsigned_type;

  // std::_Integer_class<int>::_Signed_type qqq;

  constexpr bool defint = std::default_initializable<xpm::voxel_idx>;
  constexpr bool copyable = std::copyable<xpm::voxel_idx>;
  constexpr bool movable = std::movable<xpm::voxel_idx>;

  constexpr bool check00 = std::_Integer_like<xpm::voxel_idx>;

  // constexpr bool lessthen = static_cast<xpm::voxel_idx>(-1) < static_cast<xpm::voxel_idx>(0);

  using iter_diff = std::iter_difference_t<int>;
  iter_diff lul;

  constexpr bool check1 = std::_Signed_integer_like<iter_diff>;

  constexpr bool check = std::weakly_incrementable<xpm::voxel_idx>;

  constexpr bool semireg = std::semiregular<xpm::voxel_idx>;

  // std::_Weakly_equality_comparable_with<

  // std::ranges::iota_view<int, int>{0, 20};
  

  

  for (auto rrrr : std::ranges::subrange(xpm::voxel_idx{0}, xpm::voxel_idx{20})) {
    
  }

  std::ranges::iota_view<xpm::voxel_idx, xpm::voxel_idx>{xpm::voxel_idx{0}, xpm::voxel_idx{20}};
  // std::weakly_incrementable

  // constexpr xpm::image_idx qwe3{3};
  //
  //
  //
  // constexpr auto qwewww = qwe + 5;
  //
  // constexpr auto omsodfsf = *qwewww;


  if (argc == 2 && !std::strcmp(argv[1], "-s")) {
    using namespace boost::interprocess;
    
    MPI_Init(&argc, &argv);
  
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    static constexpr auto root = 0;
    
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
    
    dpl::hypre::ls_unknown_storage lus(count, jlower);

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpi_rank == root)
    //   std::cout << "\n\nPRE_SOLVE" << std::flush;
    // MPI_Barrier(MPI_COMM_WORLD);

    // try {

    

    dpl::hypre::solve(
      lk_ref/*.get_ref()*/,
      lus.get_ref());


    std::unique_ptr<double[]> pressure_utpr;
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

    MPI_Gatherv(
      lus.data[dpl::hypre::keys::value].get(), count, MPI_DOUBLE,
      recvbuf, recvcounts, displs, MPI_DOUBLE,
      root, MPI_COMM_WORLD);

    if (mpi_rank == root) {
      shared_memory_object smo_output{open_or_create, "xpm-hypre-output", read_write};
      smo_output.truncate(lk_ref.nrows*sizeof(double));
      mapped_region region_output(smo_output, read_write);
      std::memcpy(region_output.get_address(), recvbuf, lk_ref.nrows*sizeof(double));
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