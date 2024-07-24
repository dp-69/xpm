/*
 * This file is part of Extensive Pore Modelling (xpm).
 *   | https://github.com/dp-69/xpm
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * xpm is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * xpm is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with xpm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

  // using namespace std::string_literals;

  if (argc == 2 && !std::strcmp(argv[1], "-s")) {
    MPI_Init(&argc, &argv);
    dpl::hypre::process();
    MPI_Finalize();
    return 0;
  }

  #ifdef _WIN32
    MPI_Init(&argc, &argv);
  #endif

  /*
   *
   */
  constexpr auto author_note = 
R"( * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                Extensive Pore Modelling - xpm v0.2.3                *
 *                                                                     *
 *                        Copyright (c) 2024                           *
 *   Dmytro Petrovskyy, Julien Maes, Hannah P. Menke, Kamaljit Singh   *
 *                                                                     *
 *                    https://github.com/dp-69/xpm                     *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

)";

  std::cout << author_note;

  /*
   *
   */


  dpl::mpi::exec = argv[0];
  auto cmdl = argh::parser(argc, argv);

  std::filesystem::path input;
  cmdl(1, "config.json") >> input;

  try {
    if (cmdl["G"]) { /*argc == 2 && !std::strcmp(argv[1], "-G")*/
      using json = nlohmann::json;

      auto non_gui_exec = [&input](json& j) {
        xpm::modeller modeller;
    
        modeller.init(input, j);

        modeller.prepare();
        modeller.compute_pressure();

        auto dir =
          std::filesystem::path(dpl::mpi::exec)
            .replace_filename("results")/modeller.cfg().image.path.stem();

        create_directories(dir);

        {
          nlohmann::json pp_j;
          pp_j["total_poro"] = modeller.petrophysics_summary().total_porosity;
          pp_j["total_perm"] = modeller.petrophysics_summary().perm_total;
          pp_j["macro_perm"] = modeller.petrophysics_summary().perm_macro;

          std::ofstream{dir/"petrophysics_summary.json"} << pp_j.dump(2);
        }

        if (modeller.cfg().report.invasion_percolation) {
          modeller.get_invasion_task().init();
      
          modeller.get_invasion_task().launch_primary(
            modeller.absolute_rate(),
            modeller.cfg().theta);

          std::ofstream{dir/"pc_primary.txt"} << modeller.pc_to_plain<true>();
          std::ofstream{dir/"pc_secondary.txt"} << modeller.pc_to_plain<false>();

          std::ofstream{dir/"kr_primary.txt"} << modeller.kr_to_plain<true>();
          std::ofstream{dir/"kr_secondary.txt"} << modeller.kr_to_plain<false>();

          std::cout << '\n';
        }

        std::ofstream{dir/"phi_k_kr_pc.json"} << modeller.petrophysics_json().dump(2);
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