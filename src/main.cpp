#include "xpm_widget.hpp"
#include "xpm_widget.hpp"

#include <QWidget>


int main(int argc, char* argv[])
{
  // auto j_complete = nlohmann::json::parse("[[0.018006, 9.853437e-01], [0.093010, 8.867584e-01], [0.168171, 7.531006e-01], [0.243182, 6.618532e-01], [0.318300, 5.582312e-01], [0.393701, 4.511272e-01], [0.468938, 4.274975e-01], [0.544174, 3.588589e-01], [0.620258, 2.304018e-01], [0.695897, 7.069218e-02], [0.771625, 6.008966e-02], [0.846732, 5.907953e-08], [0.923965, 2.416029e-09], [1.000000, 0.000000e+00]]");
  //
  // std::vector<dpl::vector2d> foo = j_complete;
  //
  // auto q = solve(std::span<const dpl::vector2d>{foo}, 1.0, dpl::extrapolant::flat);
  //
  // int p = 3;

  // int p = 3;
  //
  // xpm::crop(
  //   R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Est-v0m2s3_500x500x500_4p0um.raw)", 500, 128,
  //   R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Est-v0m2s3_256x256x256_4p0um.raw)", 256);


  // auto v0 = xpm::invaded_func<false>::flipper_::flip(false);
  // auto v1 = xpm::invaded_func<false>::flipper_::flip(true);
  // auto v2 = xpm::invaded_func<true>::flipper_::flip(false);
  // auto v3 = xpm::invaded_func<true>::flipper_::flip(true);

  // std::span<dpl::vector2d> s;
  // std::span<const dpl::vector2d> s_cnst;
  // std::span<const dpl::vector2d> s_cnst2;
  // s_cnst = s;

  // static const auto input = std::initializer_list<dpl::vector2d>{{1, 20}, {3, 30}, {3.0001, 40}, {7, 50}, {10, 100}};

  // std::array<int, 5> input = {int{1}, {3}, {3}, {7}, {10}};


  //  std::array<dpl::vector2d, 5> input = {dpl::vector2d{1, 20}, {3, 30}, {3.000, 40}, {7, 50}, {10, 100}};
  //
  //
  //
  //
  // auto ret = std::ranges::unique(input,
  //   {}    /*[](double l, double r) { return false; }*/,
  //   [](const dpl::vector2d& v) { return v.x(); }).begin();
  //
  // // for (auto x : foo) {
  // //   std::cout << x << '\n';
  // // }
  //
  // std::vector www (input.begin(), ret);
  //
  // int p = 3;

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
