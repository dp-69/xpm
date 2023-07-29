#include "xpm_widget.hpp"

#include <QWidget>


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  if (argc == 2 && !std::strcmp(argv[1], "-s"))
    dpl::hypre::mpi::process();
  else {
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
  }

  return MPI_Finalize();
}
