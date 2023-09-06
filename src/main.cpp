#include "xpm_widget.hpp"

#include <QWidget>


int main(int argc, char* argv[])
{
  // using BoostGraph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>;
  //
  // BoostGraph g(5);
  //
  // add_edge(1, 0, g);
  // add_edge(2, 0, g);
  //
  // auto [iter, end] = out_edges(0, g);
  // for (; iter != end; ++iter) {
  //
  //   
  //   std::cout << target(*iter, g) << ' ';
  //   
  // }
  //
  // int k = 3;

  // {
  //   using namespace std::chrono;
  //
  //   auto last = system_clock::now() - seconds{4};
  //
  //   if (auto diff = duration_cast<seconds>(system_clock::now() - last); diff < seconds{5})
  //     std::this_thread::sleep_for(seconds{5} - diff);
  // }

  // auto diff = std::chrono::system_clock::now() - 
  // (std::chrono::system_clock::now() - std::chrono::minutes{5});
  //
  //
  //
  // std::cout << duration_cast<std::chrono::minutes>(diff);
  //
  // int  p = 3;

  // {
  //   using namespace xpm;
  //
  //   displ_queue ms;
  //   ms.insert(displ_elem::macro, 1, 20.);
  //   ms.insert(displ_elem::macro, 3, 50.);
  //   ms.insert(displ_elem::macro, 2, 30.);
  //   ms.insert(displ_elem::macro, 4, 50.);
  //
  //   auto first = ms.front();
  //   ms.pop();
  //
  //   int p = 3;
  // }

  // xpm::test::DFS_CHECK();
  // xpm::test::split_join_validity_check_ET_ONLY();
  // xpm::test::split_join_validity_check_ETNTE_ONLY();
  //
  // getchar();
  // return 0;

  if (argc == 2 && !std::strcmp(argv[1], "-s")) {
    MPI_Init(&argc, &argv);
    dpl::hypre::mpi::process();
    MPI_Finalize();
    return 0;
  }
  else {
    #ifdef _WIN32
      MPI_Init(&argc, &argv);
    #endif

    dpl::hypre::mpi::mpi_exec = argv[0];

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

    /*auto result = */QApplication::exec();

    #ifdef _WIN32
      MPI_Finalize();
    #endif
  }

  return 0;
}
