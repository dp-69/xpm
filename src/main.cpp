#include "xpm_widget.hpp"
#include "xpm_widget.hpp"

#include <QWidget>


int main(int argc, char* argv[])
{
  // std::span<dpl::vector2d> s;
  // std::span<const dpl::vector2d> s_cnst;
  // std::span<const dpl::vector2d> s_cnst2;
  // s_cnst = s;

  // static const auto input = std::initializer_list<dpl::vector2d>{{1, 20}, {3, 30}, {3.0001, 40}, {7, 50}, {10, 100}};

  // constexpr std::array<dpl::vector2d, 5> input = {dpl::vector2d{1, 20}, {3, 30}, {3.0001, 40}, {7, 50}, {10, 100}};
  //
  // auto spn = std::span<const dpl::vector2d>{input};
  // auto qqq = spn.size();
  // auto val = solve(std::span<const dpl::vector2d>{input}, 0.95, dpl::extrapolant::flat);

  // auto k = solve(std::span{input}, 3);

  // xpm::idx1d_t size = 10;
  //
  // std::size_t q = 0;
  //
  //
  // for (auto qwe : std::ranges::iota_view{q, q} | std::views::filter([](std::size_t i) { return i%2; })) {
  //   std::cout << qwe << '\n';
  // }

  // std::ranges::min(pn_.node_.range(attribs::r_ins), {},
  //           [theta](double r_ins) { return props::r_cap_piston_with_films(theta, r_ins); })*0.95

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
