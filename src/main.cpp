#include "widget.hpp"

#include <QWidget>
#include <QApplication>

#include <argh.h>

int main(int argc, char* argv[])
{
  // xpm::transform(
  //   R"(C:\dev\.temp\images\176-F8-15-97_1000x500x500_29p0um.raw)",
  //   {1000, 500, 500}, 0, 
  //   R"(C:\dev\.temp\images\176-F8-15-97_500x250x250_29p0um.raw)",
  //   {500, 250, 250},
  //   [](std::uint8_t v) { return /*v < 50 ? 0 : */v; });
  //
  // getchar();
  
  // for (auto lul : std::ranges::iota_view(xpm::voxel_t{0}, 10)) {
  //   std::cout << *lul;
  // }

  // int p = 3;

  if (argc == 2 && !std::strcmp(argv[1], "-s")) {
    MPI_Init(&argc, &argv);
    dpl::hypre::process();
    MPI_Finalize();
    return 0;
  }

  #ifdef _WIN32
    MPI_Init(&argc, &argv);
  #endif

  dpl::mpi::exec = argv[0];
  auto cmdl = argh::parser(argc, argv);

  std::filesystem::path input;
  cmdl(1, "config.json") >> input;

  try {
    if (cmdl["G"]) { /*argc == 2 && !std::strcmp(argv[1], "-G")*/
      using json = nlohmann::json;

      auto non_gui_exec = [](const json& j) {
        xpm::modeller modeller;
    
        modeller.init(j);

        modeller.compute_pressure();

        if (modeller.settings().report.invasion_percolation) {
          modeller.get_invasion_task().init();
      
          modeller.get_invasion_task().launch_primary(
            modeller.absolute_rate(),
            modeller.settings().theta,
            modeller.settings().primary.pc.inverse_unique());

          auto dir =
            std::filesystem::path(dpl::mpi::exec)
              .replace_filename("results")/modeller.settings().image.path.stem();

          create_directories(dir);

          std::ofstream(dir/"pc_primary.txt") << modeller.pc_to_plain(std::true_type{});
          std::ofstream(dir/"pc_secondary.txt") << modeller.pc_to_plain(std::false_type{});

          std::ofstream(dir/"kr_primary.txt") << modeller.kr_to_plain(std::true_type{});
          std::ofstream(dir/"kr_secondary.txt") << modeller.kr_to_plain(std::false_type{});
        }

        std::cout << '\n';
        // add_copy("Copy (primary)", [this] { return modeller.pc_to_plain(std::true_type{}); });
      };

      if (auto root = json::parse(std::ifstream{input}, nullptr, true, true);
        root.is_array())
        std::ranges::for_each(root, non_gui_exec);
      else
        non_gui_exec(root);
    }
    else {
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

      QApplication app(argc, argv);

      xpm::Widget widget;
        
      widget.Init(input);
        
      // Ui::MainWindow ui;
      // ui.setupUi(&widget);

      widget.resize(1400, 1000);
      widget.show();

      /*auto result = */QApplication::exec();
    }
  }
  catch (const xpm::config_exception& e) {
    fmt::print("\n(error) {}", e.what());
  }

  #ifdef _WIN32
    MPI_Finalize();
  #endif
    
  return 0;
}




  // using namespace std::numbers;
  // using eq_tr = xpm::hydraulic_properties::equilateral_triangle_properties;
  //
  // auto theta_rec = 10*pi/180;
  // auto theta_adv = 55*pi/180; 
  //
  // auto r_ins = 0.001;
  // double drainage_ratio = 0.95;
  //
  // auto r_pld_valv = eq_tr::r_cap_piston_with_films_valvatne(theta_rec, r_ins);
  // auto pc_pld_val = 1/r_pld_valv;
  //
  // {
  //   constexpr auto G = eq_tr::area()/eq_tr::sqr(eq_tr::perimeter());
  //
  //   auto diff = r_pld_valv - eq_tr::r_cap_piston_with_films(theta_rec, r_ins);
  //
  //   auto a_films_new = eq_tr::area_corners_valv(theta_rec, r_pld_valv);
  //   auto a_films_old = eq_tr::area_of_films(theta_rec, r_pld_valv);
  //   auto a_diff = a_films_new - a_films_old;
  //
  //   auto simple_bal = eq_tr::simple_balance(theta_rec, r_pld_valv, r_ins);
  //   auto pinned_bal = eq_tr::pinned_balance(eq_tr::b_length(theta_rec, r_pld_valv), theta_adv, r_pld_valv, r_ins);
  //
  //
  //   // auto diff_snap_off = eq_tr::r_cap_collapse(theta_adv, r_ins) - eq_tr::r_cap_snap_off_valv(theta_adv, r_ins);
  //
  //   // int p = 3;
  //
  // }
  //
  // auto r_cap_min = r_pld_valv*drainage_ratio;
  // auto pc_max = 1/r_cap_min;
  //
  // auto b_rec = eq_tr::b_length(theta_rec, r_cap_min);
  //
  // auto pc_for_hinge_adv = 1/eq_tr::r_cap_hinging(b_rec, theta_adv);
  // auto pc_for_hinge_rec = 1/eq_tr::r_cap_hinging(b_rec, theta_rec);
  //
  // auto func = [=](double pc) { return eq_tr::pinned_balance(b_rec, theta_adv, 1/pc, r_ins); };
  //
  // for (auto steps = 20, i = 0; i <= steps; i++) {
  //   auto pc = pc_for_hinge_rec - i*(pc_for_hinge_rec - pc_for_hinge_adv)/steps;
  //
  //   auto theta_h = eq_tr::hinging(b_rec, 1/pc);
  //
  //   std::cout << fmt::format("pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}, hinge: {:.1f} deg, sw: {:.5f}\n",
  //     pc, 1/pc, func(pc), theta_h/pi*180, eq_tr::area_corners_valv(theta_h, 1/pc)/eq_tr::area(r_ins)); // decreases pc, and decreases balance
  // }
  //
  // // try
  // {
  //   uintmax_t max = 64;
  //
  //   auto fa = func(pc_for_hinge_adv);
  //   auto fb = func(pc_for_hinge_rec);
  //
  //   if (fa*fb < 0) {
  //     using namespace boost::math::tools;
  //
  //     auto pc = toms748_solve(
  //       func,
  //       pc_for_hinge_adv,
  //       pc_for_hinge_rec,
  //       fa,
  //       fb,
  //       eps_tolerance<double>(),
  //       max).first;
  //
  //     std::cout << fmt::format("\nSOLUTION pc: {:7.2f}, r_cap {:.3e}, iters: {}, energy: {: .3e}, hinge: {:.1f} deg\n\n",
  //       pc, 1/pc, max, func(pc), eq_tr::hinging(b_rec, 1/pc)/pi*180);
  //     
  //     // std::cout << fmt::format("\npc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}\n", result.first, 1/result.first, func(result.first)); // decreases pc, and decreases balance
  //     // std::cout << fmt::format("pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}\n", result.first/2, 2/result.first, func(result.first/2)); // decreases pc, and decreases balance
  //
  //     // auto result = boost::math::tools::bracket_and_solve_root(
  //     //   func, 
  //     //   pc_for_hinge_rec,
  //     //   2.,
  //     //   true,
  //     //   boost::math::tools::eps_tolerance<double>(), max);
  //   }
  //   else {
  //     std::cout << fmt::format("\nNO SOLUTION\n\n");
  //   }
  //
  //   {
  //     auto r_cap = eq_tr::r_cap_piston_with_films_valvatne(theta_adv, r_ins);
  //     std::cout << fmt::format("\nTRIVIAL_FILMS pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}, theta_adv: {:.1f} deg\n\n",
  //       1/r_cap, r_cap, eq_tr::simple_balance(theta_adv, r_cap, r_ins), theta_adv/pi*180);
  //   }
  //
  //   // {
  //   //   auto r_cap = eq_tr::r_cap_pison_no_films(theta_adv, r_ins);
  //   //   std::cout << fmt::format("\nTRIVIAL_NO_FILMS pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}, theta_adv: {:.1f} deg\n\n",
  //   //       1/r_cap, r_cap, eq_tr::simple_balance(theta_adv, r_cap, r_ins), theta_adv/pi*180);
  //   // }
  //
  //   {
  //     auto r_cap = eq_tr::r_cap_piston_secondary(theta_rec, r_cap_min, theta_adv, r_ins);
  //     std::cout << fmt::format("\nFINAL_FUNC pc: {:7.2f}, r_cap {:.3e}\n\n",
  //         1/r_cap, r_cap/*, eq_tr::hinging(b_rec, r_cap)/pi*180*/);
  //   }
  //      
  //
  //   getchar();
  // }



// if (pc_for_hinge_adv < 0) {
  //   std::cout << "---\n";
  //   for (auto steps = 50, i = 0; i < steps; i++) {
  //     auto pc = 0 - (i + 1)*(0 - pc_for_hinge_adv)/steps;
  //     std::cout << fmt::format("pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}\n", pc, 1/pc, func(pc)); // decreases pc, and decreases balance
  //   }
  //
  //   uintmax_t max = 64;
  //
  //   auto result = boost::math::tools::bracket_and_solve_root(
  //     func, 
  //     -0.0000001,
  //     2.,
  //     true,
  //     boost::math::tools::eps_tolerance<double>(), max);
  //
  //   std::cout << fmt::format("\npc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}\n", result.first, 1/result.first, func(result.first)); // decreases pc, and decreases balance
  //   std::cout << fmt::format("pc: {:7.2f}, r_cap {:.3e}, energy: {: .3e}\n", result.first/2, 2/result.first, func(result.first/2)); // decreases pc, and decreases balance
  //   std::cout << fmt::format("SOLUTION pc: {:7.2f}, r_cap {:.3e}, iters: {}, hinge: {:.1f} deg\n\n",
  //     result.first, 1/result.first, max, eq_tr::hinging(b_rec, 1/result.first)/pi*180);
  //
  //    getchar();
  // }
  //
  //
  //
  // getchar();









// {
  //       auto invThetaAdv = phase_properties<InvPhase>::phase_contact_angle(elem.phase1ContactAngle[1]);        
  //
  //       auto pinnedAlpha = equilateral_triangle_properties::energy_pinned_to_single(invThetaAdv, elem.drainageRatio, invTheta0).find_alpha();
  //
  //       auto pinnedToSingle = phase_properties<DefPhase>::sign*equilateral_triangle_properties::rCap(elem.abLength, pinnedAlpha);
  //
  //       if (invThetaAdv < M_PI/3.0) {
  //         if (invTheta0 <= (M_PI/3.0 - pinnedAlpha) <= invThetaAdv)
  //           queue.insert(totalIdx, pinnedToSingle, phase_distribution(InvPhase)); 
  //         else
  //           queue.insert(
  //             totalIdx, 
  //             phase_properties<DefPhase>::sign*equilateral_triangle_properties::rCap_PLD_with_wetting_films(invThetaAdv, elem.rInsEffectivePc),
  //             phase_distribution(InvPhase));
  //       }
  //       else if (2.0*M_PI/3.0 < invThetaAdv) {        
  //         auto pinnedToLayers = phase_properties<InvPhase>::sign*equilateral_triangle_properties::rCap_PLD_with_wetting_films(M_PI - invThetaAdv, elem.rInsEffectivePc);
  //         
  //         if (invasion_comparator<InvPhase>::compare(pinnedToLayers, pinnedToSingle)) {            
  //           auto layersToSingle = phase_properties<DefPhase>::sign*
  //             equilateral_triangle_properties::rCap(elem.abLength, 
  //               equilateral_triangle_properties::energy_layers_to_single(invThetaAdv).find_alpha());
  //
  //           queue.insert(totalIdx, pinnedToLayers, phase_distribution(InvPhase, pore_occupancy_type::blc));
  //           queue.insert(totalIdx, layersToSingle, phase_distribution(InvPhase));
  //
  //         }            
  //         else
  //           queue.insert(totalIdx, pinnedToSingle, phase_distribution(InvPhase));
  //       }
  //       else
  //         queue.insert(totalIdx, pinnedToSingle, phase_distribution(InvPhase));
  //     }










