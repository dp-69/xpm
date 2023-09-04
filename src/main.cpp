#include "xpm_widget.hpp"

#include <QWidget>


int main(int argc, char* argv[])
{
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
