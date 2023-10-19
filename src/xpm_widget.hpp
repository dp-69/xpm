#pragma once


#include "functions.h"
#include "displ_queue.hpp"
#include "invasion_task.hpp"


#include <dpl/units.hpp>
#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/hypre/mpi_module.hpp>
#include <dpl/qt/layout.hpp>
#include <dpl/qt/property_editor/PropertyItemsBase.hpp>
#include <dpl/qt/property_editor/QPropertyTreeView.hpp>
#include <dpl/vtk/ImageDataGlyphMapper.hpp>
#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>
  
// #include <QWidget>
#include <QMainWindow>
#include <QSplitter>
#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>
#include <QLogValueAxis>
#include <QClipboard>
#include <QLegendMarker>

#if (VTK_MAJOR_VERSION == 8)
  #include <QVTKOpenGLWidget.h>
  #include <QSurfaceFormat>
#elif (VTK_MAJOR_VERSION == 9)
  #include <QVTKOpenGLNativeWidget.h>
#endif


#include <vtkAssembly.h>
#include <vtkCellData.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGlyph3DMapper.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <algorithm>
#include <future>
#include <unordered_set>

#include <nlohmann/json.hpp>

// #undef LoadImage

namespace xpm
{
  #if (VTK_MAJOR_VERSION == 8)
    using QVTKWidgetRef = QVTKOpenGLWidget;
  #elif (VTK_MAJOR_VERSION == 9)
    using QVTKWidgetRef = QVTKOpenGLNativeWidget;
  #endif


  class XPMWidget : public QMainWindow
  {
  Q_OBJECT

    using QPropertyTreeView = dpl::qt::property_editor::QPropertyTreeView;
    
    QPropertyTreeView* tree_view_ = nullptr; 
    
    vtkNew<vtkRenderer> renderer_; 
    vtkNew<vtkGenericOpenGLRenderWindow> render_window_;
    
    // vtkRenderWindowInteractor* interactor_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;

    vtkNew<vtkLookupTable> lut_pressure_;
    vtkNew<vtkLookupTable> lut_continuous_;
    // vtkNew<vtkLookupTable> lut_velem_;
    vtkNew<vtkLookupTable> lut_pore_solid_micro_;
    vtkNew<vtkLookupTable> lut_node_throat_;

    vtkFloatArray* macro_colors;
    vtkFloatArray* throat_colors;


    dpl::vtk::ImageDataGlyphMapper<idx1d_t> img_glyph_mapper_;






    pore_network pn_;
    image_data img_;
    pore_network_image pni_{pn_, img_};
    // occupancy_arrays occupancy_arrays_;

    std::array<double, 6> bounds_ = {0, 100, 0, 100, 0, 100};

    startup_settings settings_;

    invasion_task invasion_task_{pni_, settings_};


    std::string status_ = "<nothing>";
    dpl::qt::property_editor::PropertyItem* status_property_item_;

    QChartView* pc_chart_view_;
    QChartView* kr_chart_view_;
    QLineSeries* line_series_;
    QLineSeries* kr0_series_;
    QLineSeries* kr1_series_;
    double absolute_rate_;


    QValueAxis* kr_axis_y_;
    QLogValueAxis* pc_axis_y_;
    
    void UpdateStatus(std::string text) {
      status_ = std::move(text);
      auto idx = tree_view_->model()->index(status_property_item_, 1);
      tree_view_->model()->dataChanged(idx, idx);
    }

    std::future<void> consumer_future_;
    





    


    
    
    

    static void InitLutNodeThroat(vtkLookupTable* lut) {
      lut->IndexedLookupOn();
      lut->SetNumberOfTableValues(2);
      lut->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut->SetAnnotation(vtkVariant(0), "Node");
      lut->SetTableValue(1, 0.784314, 0.752624, 0.467419);
      lut->SetAnnotation(vtkVariant(1), "Throat");
    }

    static void InitLutPoreSolidMicro(vtkLookupTable* lut) {
      lut->IndexedLookupOn();  
      lut->SetNumberOfTableValues(3);
      lut->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut->SetAnnotation(vtkVariant(0), "Pore");
      lut->SetTableValue(1, 0.5, 0.5, 0.5);
      lut->SetAnnotation(vtkVariant(1), "Solid");
      lut->SetTableValue(2, 0.8, 0.8, 0.8);
      lut->SetAnnotation(vtkVariant(2), "Microporous");
    }

    static void InitLutVelems(vtkLookupTable* lut, std::int32_t max) {
      auto count = max - (-2)/**min*/ + 1;

      lut->SetNumberOfTableValues(count);
      
      for (int32_t i = 0; i < count; ++i) {
        auto coef = static_cast<double>(i*45%count)/count;
        auto color = QColor::fromHsl(coef*255, 175, 122);  // NOLINT(cppcoreguidelines-narrowing-conversions)
        lut->SetTableValue(i, color.redF(), color.greenF(), color.blueF());  // NOLINT(clang-diagnostic-double-promotion)
      }

      lut->SetTableValue(0, 0.7, 0.7, 0.7);     // -2 velem value
      lut->SetTableValue(1, 0.7, 0.7, 0.7);     // -1 velem value
      lut->SetTableRange(-2/**min*/, max);
    }


    

    void InitGUI() {
      qvtk_widget_ = new QVTKWidgetRef;

      #if (VTK_MAJOR_VERSION == 8)
        qvtk_widget_->SetRenderWindow(render_window_);
      #elif (VTK_MAJOR_VERSION == 9)
        qvtk_widget_->setRenderWindow(render_window_);
      #endif
      
      render_window_->AddRenderer(renderer_);
      // interactor_ = render_window_->GetInteractor();

      renderer_->SetBackground(dpl::vector3d{1});

      {
        tree_view_ = new QPropertyTreeView;

        auto* model = tree_view_->model();
        
        model->AddItem(
          dpl::qt::property_editor::ItemFunctor<bool>{
            "Edges",
            [this] {
              return static_cast<bool>(std::get<0>(img_glyph_mapper_.faces_).GetActor()->GetProperty()->GetEdgeVisibility());
            },
            [this](bool v) {
              dpl::sfor<6>([this, v](auto i) {
                static_cast<dpl::vtk::GlyphMapperFace<idx1d_t>&>(std::get<i>(img_glyph_mapper_.faces_))
                  .GetActor()->GetProperty()->SetEdgeVisibility(v);
                render_window_->Render();
              });
            }
          });

        model->AddItem(
          dpl::qt::property_editor::ItemFunctor<bool>{
            "Microporosity",
            [this] {
              return static_cast<bool>(std::get<0>(img_glyph_mapper_.faces_).GetActor()->GetVisibility());
            },
            [this](bool v) {
              dpl::sfor<6>([this, v](auto i) {
                static_cast<dpl::vtk::GlyphMapperFace<idx1d_t>&>(std::get<i>(img_glyph_mapper_.faces_))
                  .GetActor()->SetVisibility(v);
                render_window_->Render();
              });
            }
          });

        status_property_item_ = model->AddItem(
          dpl::qt::property_editor::ItemFunctor<std::string>{"Status", [this] { return status_; }});

        auto* plot_tabs = new QTabWidget;

        {
          pc_chart_view_ = new QChartView;
          pc_chart_view_->setRenderHint(QPainter::Antialiasing);
          pc_chart_view_->setBackgroundBrush(Qt::GlobalColor::white);

          {
            auto menu = std::make_unique<QMenu>();
            
            auto action = [this](auto format, auto sep, auto begin, auto end) {
              return [=, this] {
                auto view =
                  std::views::transform(invasion_task_.pc_curve(),
                    [format](const dpl::vector2d& p) { return fmt::format(format, p.x(), p.y()); })
                  | std::views::reverse;

                std::stringstream ss;
                ss << begin << *std::begin(view);
                for (const auto& s : view | std::views::drop(1))
                  ss << sep << s;
                ss << end << '\n';

                QApplication::clipboard()->setText(QString::fromStdString(ss.str()));
              };
            };

            using fmt_string = fmt::format_string<const double&, const double&>;

            connect(menu->addAction("Copy"), &QAction::triggered,
              action(fmt_string{"{:.6f}\t{:.6e}"}, "\n", "", ""));

            connect(menu->addAction("Copy JSON"), &QAction::triggered,
              action(fmt_string{"[{:.6f}, {:.6e}]"}, ", ", "[", "]"));

            pc_chart_view_->setContextMenuPolicy(Qt::ContextMenuPolicy::CustomContextMenu);
            connect(pc_chart_view_, &QTreeView::customContextMenuRequested,
              [this, m = std::move(menu)](const QPoint& pos) { m->popup(pc_chart_view_->viewport()->mapToGlobal(pos)); });
          }
        
          auto* chart = new QChart;
          pc_chart_view_->setChart(chart);

          auto* axis_x = new QValueAxis;//static_cast<QValueAxis*>(chart_->axes(Qt::Horizontal)[0]);
          axis_x->setLabelFormat("%.2f");
          axis_x->setTitleText("Water saturation");
          axis_x->setRange(0.0, 1);
          chart->addAxis(axis_x, Qt::AlignBottom);

          pc_axis_y_ = new QLogValueAxis;//static_cast<QLogValueAxis*>(chart_->axes(Qt::Vertical)[0]);
          pc_axis_y_->setLabelFormat("%.0e");
          pc_axis_y_->setTitleText("Capillary pressure, Pa");
          chart->addAxis(pc_axis_y_, Qt::AlignLeft);

          line_series_ = new QLineSeries;
          line_series_->setName("total");
          line_series_->clear();
          line_series_->setPointsVisible(true);

          chart->addSeries(line_series_);
          line_series_->attachAxis(axis_x);
          line_series_->attachAxis(pc_axis_y_);

          plot_tabs->addTab(pc_chart_view_, "Pc");
        }

        {
          kr_chart_view_ = new QChartView;
          kr_chart_view_->setRenderHint(QPainter::Antialiasing);
          kr_chart_view_->setBackgroundBrush(Qt::GlobalColor::white);

          {
            auto menu = std::make_unique<QMenu>();
            
            auto text = [this](auto format, auto sep, auto begin, auto end) {
              auto view =
                std::views::transform(invasion_task_.kr_curves(),
                  [format](const dpl::vector3d& p) { return fmt::format(format, p.x(), p.y(), p.z()); })
                | std::views::reverse;

              std::stringstream ss;
              ss << begin << *std::begin(view);
              for (const auto& s : view | std::views::drop(1))
                ss << sep << s;
              ss << end;

              return ss.str();
            };

            using fmt_string = fmt::format_string<const double&, const double&, const double&>;

            connect(menu->addAction("Copy"), &QAction::triggered,
              [text] {
                QApplication::clipboard()->setText(QString::fromStdString(
                  text(fmt_string{"{:.6f}\t{:.6e}\t{:.6e}"}, "\n", "", "\n")));
              });

            connect(menu->addAction("Copy JSON"), &QAction::triggered,
              [text] {
                QApplication::clipboard()->setText(QString::fromStdString(
                  fmt::format(
                    "[\n"
                    "  {},\n"
                    "  {}\n"
                    "]\n",
                    text(fmt_string{"[{0:.6f}, {1:.6e}]"}, ", ", "[", "]"),
                    text(fmt_string{"[{0:.6f}, {2:.6e}]"}, ", ", "[", "]"))
                ));
              });

            kr_chart_view_->setContextMenuPolicy(Qt::ContextMenuPolicy::CustomContextMenu);
            connect(kr_chart_view_, &QTreeView::customContextMenuRequested,
              [this, m = std::move(menu)](const QPoint& pos) { m->popup(kr_chart_view_->viewport()->mapToGlobal(pos)); });
          }

          auto* chart = new QChart;
          kr_chart_view_->setChart(chart);

          auto* axis_x = new QValueAxis;
          axis_x->setLabelFormat("%.2f");
          axis_x->setTitleText("Water saturation");
          axis_x->setRange(0, 1);
          chart->addAxis(axis_x, Qt::AlignBottom);

          kr_axis_y_ = new QValueAxis;
          kr_axis_y_->setLabelFormat("%.2f");
          kr_axis_y_->setTitleText("Relative permeability");
          kr_axis_y_->setRange(0, 1);
          chart->addAxis(kr_axis_y_, Qt::AlignLeft);

          kr0_series_ = new QLineSeries;
          kr0_series_->setPointsVisible(true);
          kr0_series_->setName("water");

          kr1_series_ = new QLineSeries;
          kr1_series_->setPointsVisible(true);
          kr1_series_->setName("phase 1");
          // line_series->setPointsVisible(true);

          // for (auto size = 50, i = 0; i <= size; ++i) {
          //   auto So = 1.*i/size;
          //   auto Sor = 0.0;
          //   auto Lambda = 2.0;
          //   auto kro = std::pow((So - Sor)/(1 - Sor), (2 + 3*Lambda)/Lambda);
          //   auto krg = std::pow((1 - So)/(1 - Sor), 2)*(1 - std::pow((So - Sor)/(1 - Sor), (2 + Lambda)/Lambda));
          //
          //   kro_series->append(So, kro);
          //   krg_series->append(So, krg);
          // }
          
          // pow((So - Sor)/(1 - Sor), (2 + 3*Lambda)/Lambda)

          // line_series->append(0.1, 0.1);
          // line_series->append(0.5, 0.8);

          chart->addSeries(kr0_series_);
          kr0_series_->attachAxis(axis_x);
          kr0_series_->attachAxis(kr_axis_y_);

          chart->addSeries(kr1_series_);
          kr1_series_->attachAxis(axis_x);
          kr1_series_->attachAxis(kr_axis_y_);

          plot_tabs->addTab(kr_chart_view_, "kr");
        }


        

        using namespace dpl::qt::layout;

        setCentralWidget(
          splitter{dir::horizontal}
          << (
            splitter{dir::vertical}
            << tree_view_ << stretch{1}
            << plot_tabs/*chart_view_*/ << stretch{0}) << stretch{0}
          << qvtk_widget_ << stretch{1}
        );
      }
    }



    auto ProcessImage(const std::filesystem::path& path) {
      namespace fs = std::filesystem;

      auto filename = path.filename();

      if (fs::path copy_path = "pnextract"/filename; absolute(copy_path) != absolute(path))
        copy(path, copy_path, fs::copy_options::update_existing);

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      current_path(prev/"pnextract");
      auto network_dir = path.stem();

      if (std::ranges::all_of(files, [&](std::string_view file) { return exists(network_dir/file); }))
        std::cout << "using cached network\n\n";
      else {
        std::cout << "=========== pnextract's network extraction begin ===========\n";

        std::system( // NOLINT(concurrency-mt-unsafe)
          fmt::format("{} {}", fs::current_path()/"pnextract", filename).c_str());
      
        create_directory(network_dir);
        for (fs::path file : files)
          rename(file, network_dir/file);

        remove(fs::path{"_VElems.mhd"});

        std::cout << "=========== pnextract's network extraction end ===========\n\n";
      }

      network_dir = absolute(network_dir);

      current_path(prev);

      return network_dir;
    }


  public:
    void Init() {
      std::locale::global(std::locale("en_US.UTF-8"));

      if (std::filesystem::exists("config.json"))
        settings_.load(nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true));

      InitGUI();

      ComputePressure();

      renderer_->ResetCamera(bounds_.data());
      renderer_->GetActiveCamera()->Zoom(0.70);
      // renderer_->ResetCameraClippingRange();
      
      tidy_axes_.Init(renderer_.Get());

      connect(qvtk_widget_, &QVTKWidgetRef::resized, this, [this] { tidy_axes_.RefreshAxes(); });

      // tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // tidy_axes_.SetFormat(".2e");
      tidy_axes_.Build(bounds_.data());

      // {
      //   net_idx_t connected = pni_.connected_total_count();
      //
      //   occupancy_arrays_.macro.resize(*connected);
      //   for (net_idx_t i{0}; i < connected; ++i)
      //     occupancy_arrays_.macro[i] = 0;
      // }

      LaunchInvasion();
    }

    void LaunchInvasion() {
      {
        auto* chart = pc_chart_view_->chart();
        auto* darcy_pc_series = new QLineSeries;
        chart->addSeries(darcy_pc_series);

        for (auto& p : settings_.darcy_pc)
          darcy_pc_series->append(p.x(), p.y());
        darcy_pc_series->setName("microporous");
        darcy_pc_series->attachAxis(chart->axes(Qt::Horizontal)[0]);
        darcy_pc_series->attachAxis(pc_axis_y_);
        chart->legend()->markers(darcy_pc_series)[0]->setVisible(false);
        darcy_pc_series->setPen(QPen{Qt::gray, 1, Qt::DashLine});
      }

      {
        auto* chart = kr_chart_view_->chart();

        auto* darcy_kr0_series = new QLineSeries;
        chart->addSeries(darcy_kr0_series);

        for (auto& p : settings_.darcy_kr0)
          darcy_kr0_series->append(p.x(), p.y());
        darcy_kr0_series->attachAxis(chart->axes(Qt::Horizontal)[0]);
        darcy_kr0_series->attachAxis(kr_axis_y_);
        chart->legend()->markers(darcy_kr0_series)[0]->setVisible(false);
        darcy_kr0_series->setPen(QPen{Qt::gray, 1, Qt::DashLine});

        auto* darcy_kr1_series = new QLineSeries;
        chart->addSeries(darcy_kr1_series);

        for (auto& p : settings_.darcy_kr1)
          darcy_kr1_series->append(p.x(), p.y());
        darcy_kr1_series->attachAxis(chart->axes(Qt::Horizontal)[0]);
        darcy_kr1_series->attachAxis(kr_axis_y_);
        chart->legend()->markers(darcy_kr1_series)[0]->setVisible(false);
        darcy_kr1_series->setPen(QPen{Qt::gray, 1, Qt::DashLine});
      }


      using eq_tr = hydraulic_properties::equilateral_triangle_properties;

      {
        using namespace std;

        auto min_r_cap_throat = numeric_limits<double>::max();

        for (size_t i{0}; i < pn_.throat_count(); ++i)
          if (auto [l, r] = pn_.throat_[attribs::adj][i]; pn_.inner_node(r) && pni_.connected(l))
            min_r_cap_throat = min(min_r_cap_throat, pn_.throat_[attribs::r_ins][i]);

        pc_axis_y_->setRange(
          // 1e4/*axis_y->min()*/,
          pow(10., floor(log10(1./eq_tr::r_cap_piston_with_films(0, ranges::max(pn_.node_.span(attribs::r_ins)))))),
          pow(10., ceil(log10(max(
            1/(0.7*eq_tr::r_cap_piston_with_films(0, min_r_cap_throat)),
            settings_.darcy_pc.empty() ? 0 : settings_.darcy_pc.front().y()))))*1.01
        );
      }

      invasion_task_.init();
      invasion_task_.generate_graph();

      consumer_future_ = std::async(std::launch::async, [this] {
        auto start = std::chrono::system_clock::now();
        double theta = 0*std::numbers::pi/180;

        using namespace std::ranges::views;
        auto query = settings_.darcy_pc | reverse | transform([](const dpl::vector2d& p) { return dpl::vector2d{p.y(), p.x()}; });
        std::vector<dpl::vector2d> pc_inv{query.begin(), query.end()};
        pc_inv.resize(std::ranges::unique(pc_inv, {}, [](const dpl::vector2d& p) { return p.x(); }).begin() - pc_inv.begin());

        using pc_sw_span = std::span<const dpl::vector2d>;




        auto invasion_future = std::async(std::launch::async, &invasion_task::launch_primary, &invasion_task_,
          absolute_rate_,
          theta,
          pc_sw_span{pc_inv});

        auto last_progress_idx = std::numeric_limits<idx1d_t>::max();

        auto update = [this, theta, start, &pc_inv, &last_progress_idx] {
          if (last_progress_idx == invasion_task_.progress_idx())
            return;

          last_progress_idx = invasion_task_.progress_idx();

          auto& state = invasion_task_.state();
          auto darcy_saturation = 1.0 - (pc_inv.empty() ? 0.0 : solve(pc_sw_span{pc_inv}, 1/state.r_cap_global, dpl::extrapolant::flat));
          
          auto map_satur = [](double x) { return x/2. + 0.25; };

          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
          
            for (vtkIdType i = 0; auto idx1d : face.GetIndices()) {
              if (auto v_idx = voxel_idx_t{idx1d};
                pni_.connected(v_idx) && state.config(pni_.net(v_idx)).phase() == phase_config::phase1())
                face.GetColorArray()->SetTypedComponent(i, 0, map_satur(darcy_saturation));
              else
                face.GetColorArray()->SetTypedComponent(i, 0, 0.0);

              ++i;
            }
          });

          using namespace attribs;

          for (macro_idx_t i{0}; i < pn_.node_count(); ++i)
            if (pni_.connected(i) && state.config(pni_.net(i)).phase() == phase_config::phase1())
              macro_colors->SetTypedComponent(*i, 0,
                map_satur(1.0 - eq_tr::area_corners(theta, state.r_cap(pni_.net(i)))/eq_tr::area(r_ins(pn_, i))));
            else
              macro_colors->SetTypedComponent(*i, 0, 0.0);
          
          {
            vtkIdType t_inner_idx = 0;

            for (std::size_t i{0}; i < pn_.throat_count(); ++i) {
              if (auto [l, r] = adj(pn_, i); pn_.inner_node(r)) {
                if (pni_.connected(l) && state.config(i).phase() == phase_config::phase1())
                  throat_colors->SetTypedComponent(t_inner_idx, 0,
                    map_satur(1.0 - eq_tr::area_corners(theta, state.r_cap(i))/eq_tr::area(r_ins(pn_, i))));
                else
                  throat_colors->SetTypedComponent(t_inner_idx, 0, 0.0);
          
                ++t_inner_idx;
              }
            }
          }
          
          QMetaObject::invokeMethod(this,
            [this, start] {
              using namespace std::chrono;

              UpdateStatus(invasion_task_.finished()
                ? fmt::format("done {}s", duration_cast<seconds>(system_clock::now() - start).count())
                : fmt::format("{:.1f} %", 100.*invasion_task_.progress_idx()/(*pni_.connected_count())));

              line_series_->clear();
              kr0_series_->clear();
              kr1_series_->clear();

              for (auto p : invasion_task_.pc_curve())
                line_series_->append(p.x(), p.y());

              for (auto [sw, kro, krg] : invasion_task_.kr_curves()) {
                kr0_series_->append(sw, kro);
                kr1_series_->append(sw, krg);                
              }

              dpl::sfor<6>([&](auto face_idx) {
                std::get<face_idx>(img_glyph_mapper_.faces_).GetColorArray()->Modified();
              });

              throat_colors->Modified();
              macro_colors->Modified();
          
              render_window_->Render();
            });
        };


        while (invasion_future.wait_for(std::chrono::milliseconds{1500}) != std::future_status::ready)
          update();
        update();
      });
    }

    void ComputePressure() {
      using clock = std::chrono::high_resolution_clock;
      using seconds = std::chrono::seconds;

      std::cout << fmt::format("image path: {}\n\n", settings_.image.path);

      auto pnm_path = ProcessImage(settings_.image.path)/"";

      auto begin_init_time = clock::now();

      //------------------------------

      dpl::vector3i processors{1};

      if (settings_.solver.decomposition)
        processors = *settings_.solver.decomposition;
      else {
        if (auto proc_count = std::thread::hardware_concurrency();
          proc_count == 12)
          processors = {2, 2, 3};
        else if (proc_count == 24)
          processors = {4, 3, 2};
        else if (proc_count == 32)
          processors = {4, 4, 2};
        else if (proc_count == 48)
          processors = {4, 4, 3};
      }

      auto cache_path = fmt::format("cache/{}-pressure-{:.2f}mD.bin",
        settings_.image.path.stem(), settings_.darcy_perm/presets::darcy_to_m2*1e3);

      InitLutNodeThroat(lut_node_throat_);
      InitLutPoreSolidMicro(lut_pore_solid_micro_);

      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      lut_pressure_->SetTableRange(0, 1);
      // lut_pressure_->SetNanColor(0.584314, 0.552624, 0.267419, 1);

      pn_.read_from_text_file(pnm_path);

      #ifdef XPM_DEBUG_OUTPUT
        std::cout
          << fmt::format("network\n  nodes: {:L}\n  throats: {:L}\n", pn_.node_count(), pn_.throat_count())
          << (pn_.eval_inlet_outlet_connectivity() ? "(connected)" : "(disconected)") << '\n';
      #endif

      #ifdef _WIN32
        pn_.connectivity_flow_summary(settings_.solver.tolerance, settings_.solver.max_iterations);
      #else
        pn_.connectivity_flow_summary_MPI(startup.solver.tolerance, startup.solver.max_iterations);
      #endif

      std::cout << '\n';

      {
        img_.read_image(settings_.image.path, settings_.image.phases);
        img_.set_dim(settings_.loaded ? settings_.image.size : std::round(std::cbrt(*img_.size())));

        img_.read_icl_velems(pnm_path);

        // auto mapped_range =
        //   std::span(img_.velem.data(), img_.size) 
        // | std::views::transform([](voxel_property::velem_t x) { return *x; });

        // InitLutVelems(lut_velem_, *std::ranges::max_element(mapped_range)); // TODO max is not valid, should check value validity

        img_.eval_microporous();
      }
      
      std::cout << "connectivity (isolated components)...";

      pni_.evaluate_isolated();

      std::cout << fmt::format(" done\n  macro: {:L}\n  voxel: {:L}\n\n",
        pni_.connected_macro_count(),
        pni_.connected_count() - pni_.connected_macro_count());

      using namespace presets;

      dpl::strong_vector<net_tag, double> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;


      if (settings_.use_cache && std::filesystem::exists(cache_path)) {
        std::cout << "using cached pressure\n\n";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(HYPRE_Complex);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(net_idx_t(nrows));
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));
      }
      else {
        std::cout << "decomposition...";

        
        auto [_, mapping] = pni_.generate_mapping(processors);

        // for (auto i = 0; i < processors.prod(); ++i)
        //   std::cout << std::format("\nblock {}, rows {}--{}, size {}",
        //     i, decomposed.blocks[i].lower, decomposed.blocks[i].upper, decomposed.blocks[i].upper - decomposed.blocks[i].lower + 1);

        std::cout
          << " done\n\n"
          << "input matrix build...";


        idx1d_t nrows = *pni_.connected_count();
        auto [nvalues, input] = pni_.generate_pressure_input(nrows, mapping.forward, single_phase_conductance{&pn_, settings_.darcy_perm});


        std::cout << " done\n\n";

        std::cout << 
          fmt::format("save input matrix [{} MB]...", (
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)))/1024/1024);

        dpl::hypre::mpi::save(input, nrows, nvalues, mapping.block_rows, settings_.solver.tolerance, settings_.solver.max_iterations);

        std::cout
          << " done\n\n"
          << "decomposition: (" << processors << fmt::format(") = {} procs\n\n", processors.prod())
          << "hypre MPI solve...";
        
        auto start = clock::now();
        
        std::system(  // NOLINT(concurrency-mt-unsafe)
          fmt::format("mpiexec -np {} \"{}\" -s", processors.prod(), dpl::hypre::mpi::mpi_exec).c_str()); 
        
        auto stop = clock::now();
        
        cout << " done " << duration_cast<std::chrono::seconds>(stop - start).count() << "s\n\n";

        {
          std::unique_ptr<HYPRE_Complex[]> decomposed_pressure;
          std::tie(decomposed_pressure, residual, iters) = dpl::hypre::mpi::load_values(nrows);

          pressure.resize(net_idx_t(nrows));

          for (HYPRE_BigInt i = 0; i < nrows; ++i)
            pressure[mapping.backward[i]] = decomposed_pressure[i];
        }

        std::cout << "pressure solved\n\n";

        if (settings_.save_cache) {
          std::filesystem::create_directory("cache");
          std::ofstream cache_stream(cache_path, std::ofstream::binary);
          cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)**pni_.connected_count());
          std::cout << "pressure cached\n\n";
        }
      }

      {
        
        auto [inlet, outlet] = pni_.flow_rates(pressure, single_phase_conductance{&pn_, settings_.darcy_perm});
          
        std::cout << fmt::format(
          "microprs perm: {:.6f} mD\n"
          "  inlet perm: {:.6f} mD\n"
          "  outlet perm: {:.6f} mD\n"
          "  residual: {:.4g}, iterations: {}\n\n",
          settings_.darcy_perm/darcy_to_m2*1000,
          inlet/pn_.physical_size.x()/darcy_to_m2*1000,
          outlet/pn_.physical_size.x()/darcy_to_m2*1000,
          residual, iters);

        absolute_rate_ = inlet;
      }
      

      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      {
        vtkSmartPointer<vtkActor> actor;

        std::tie(actor, macro_colors) = CreateNodeActor(pn_, lut_pressure_, 
          [&](macro_idx_t i) {
            return pni_.connected(i) ? /*0*/ pressure[pni_.net(i)] : std::numeric_limits<double>::quiet_NaN();
          });

        assembly->AddPart(actor);
        
        std::tie(actor, throat_colors) = CreateThroatActor(pn_, lut_pressure_, [&](std::size_t i) {
          auto [l, r] = pn_.throat_[attribs::adj][i];

          return pni_.connected(l)
            ? /*0*/ (pressure[pni_.net(l)] + pressure[pni_.net(r)])/2
            : std::numeric_limits<double>::quiet_NaN();
        });

        assembly->AddPart(actor);
      }

      renderer_->AddActor(assembly);

      // std::ranges::fill(pressure, 0);

      {
        auto scale_factor = /*1.0*/pn_.physical_size.x()/img_.dim().x(); // needed for vtk 8.2 floating point arithmetics
            
        img_glyph_mapper_.Init(scale_factor);
        {
          std::vector<bool> filter(*img_.size());
          for (voxel_idx_t i{0}; i < img_.size(); ++i)
            filter[*i] = img_.phase[i] == microporous;

          cout << "3D faces...";

          auto t0 = clock::now();

          img_glyph_mapper_.Populate(img_.dim(), pn_.physical_size/img_.dim(),
            [&](idx1d_t idx1d) { return filter[idx1d]; });

          std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
              
          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);

            idx1d_t i = 0;
            for (auto idx1d : face.GetIndices())
              face.GetColorArray()->SetTypedComponent(i++, 0, 
                pni_.connected(voxel_idx_t{idx1d})
                  ? /*0*/ pressure[pni_.net(voxel_idx_t{idx1d})]
                  : std::numeric_limits<double>::quiet_NaN()
              );

            auto* glyphs = face.GetGlyphMapper();
            auto* actor = face.GetActor();

            // glyphs->SetLookupTable(lut_velem_);
            glyphs->SetLookupTable(lut_pressure_);
            // glyphs->SetLookupTable(lut_image_phase_);
              
            glyphs->SetColorModeToMapScalars();
            glyphs->UseLookupTableScalarRangeOn();
            // glyphs->SetScalarModeToUsePointData();
            // glyphs->SetScalarModeToUseCellData();
              
            actor->SetMapper(glyphs);

            actor->GetProperty()->SetEdgeVisibility(false);
            actor->GetProperty()->SetEdgeColor(dpl::vector3d{0.25});
                
            actor->GetProperty()->SetAmbient(0.5);
            actor->GetProperty()->SetDiffuse(0.4);
            actor->GetProperty()->BackfaceCullingOn();
              
            renderer_->AddActor(actor);
          });
        }
      }


      
      tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // // tidy_axes_.SetFormat(".2e");

      bounds_ = {
        0., pn_.physical_size.x(),
        0., pn_.physical_size.y(),
        0., pn_.physical_size.z()};
    }
  };
}