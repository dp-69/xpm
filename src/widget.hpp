#pragma once

#include "modeller.hpp"


#include "functions.h"

#include <dpl/units.hpp>
#include <dpl/hypre/InputDeprec.hpp>
#include <dpl/hypre/mpi_module.hpp>
#include <dpl/qt/layout.hpp>
#include <dpl/qt/property_editor/PropertyItemsBase.hpp>
#include <dpl/qt/property_editor/QPropertyTreeView.hpp>
#include <dpl/vtk/ImageDataGlyphMapper.hpp>
#include <dpl/vtk/TidyAxes.hpp>
#include <dpl/vtk/Utils.hpp>

#include <QApplication>
#include <QChartView>
#include <QClipboard>
#include <QLegendMarker>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QMainWindow>
#include <QMenu>
#include <QTreeView>
#include <QValueAxis>

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

#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>

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


  class Widget : public QMainWindow
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

    vtkNew<vtkAssembly> macro_network_;
    vtkFloatArray* macro_colors;
    vtkFloatArray* throat_colors;


    dpl::vtk::ImageDataGlyphMapper<idx1d_t> img_glyph_mapper_;





    modeller model_;

    auto& pni() { return model_.pni(); }
    auto& pn() { return pni().pn(); }
    auto& img() { return pni().img(); }

    // pore_network pn_;
    // image_data img_;
    // pore_network_image pni_{pn_, img_};

    // runtime_settings settings_;


    std::array<double, 6> bounds_ = {0, 100, 0, 100, 0, 100};

    

    std::string status_ = "<nothing>";
    dpl::qt::property_editor::PropertyItem* status_property_item_;

    QChartView* pc_chart_view_;
    QChartView* kr_chart_view_;

    struct curves_series
    {
      QLineSeries* pc;
      QLineSeries* kr0;
      QLineSeries* kr1;
    };

    curves_series primary_;
    curves_series secondary_;

    

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
      {
        InitLutNodeThroat(lut_node_throat_);
        InitLutPoreSolidMicro(lut_pore_solid_micro_);

        dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
        dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

        lut_pressure_->SetTableRange(0, 1);
        // lut_pressure_->SetNanColor(0.584314, 0.552624, 0.267419, 1);
      }

      // this->setWindowTitle(QString::fromStdString(fmt::format("xpm - {}", settings_.image.path.filename().string())));

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

                // static_cast<dpl::vtk::GlyphMapperFace<idx1d_t>&>(std::get<i>(img_glyph_mapper_.faces_))
                //   .GetActor()->GetProperty()->SetEdgeColor(0.8, 0.8, 0.8);
                render_window_->Render();
              });
            }
          });

        model->AddItem(
          dpl::qt::property_editor::ItemFunctor<bool>{
            "Macro network",
            [this] {
              return static_cast<bool>(macro_network_->GetVisibility());
            },
            [this](bool v) {
              macro_network_->SetVisibility(v);
              render_window_->Render();
            }
          });

        model->AddItem(
          dpl::qt::property_editor::ItemFunctor<bool>{
            "Darcy nodes",
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


        tree_view_->resizeColumnToContents(0);
        // tree_view_->setColumnWidth(0, 120);


        auto* plot_tabs = new QTabWidget;

        {
          pc_chart_view_ = new QChartView;
          pc_chart_view_->setRenderHint(QPainter::Antialiasing);
          pc_chart_view_->setBackgroundBrush(Qt::GlobalColor::white);

          {
            auto menu = std::make_unique<QMenu>();

            auto add = [&](auto caption, auto text_map) {
              connect(menu->addAction(caption), &QAction::triggered, [=] {
                QApplication::clipboard()->setText(QString::fromStdString(text_map()));
              });
            };

            add("Copy (primary)", [this] { return model_.pc_to_plain(std::true_type{}); });
            add("Copy (secondary)", [this] { return model_.pc_to_plain(std::false_type{}); });

            menu->addSeparator();

            add("Copy JSON (primary)", [this] { return model_.pc_to_json(std::true_type{}); });
            add("Copy JSON (secondary)", [this] { return model_.pc_to_json(std::false_type{}); });

            pc_chart_view_->setContextMenuPolicy(Qt::ContextMenuPolicy::CustomContextMenu);
            connect(pc_chart_view_, &QTreeView::customContextMenuRequested,
              [this, m = std::move(menu)](const QPoint& pos) { m->popup(pc_chart_view_->viewport()->mapToGlobal(pos)); });
          }
        
          auto* chart = new QChart;
          pc_chart_view_->setChart(chart);

          auto* axis_x = new QValueAxis;//static_cast<QValueAxis*>(chart_->axes(Qt::Horizontal)[0]);
          axis_x->setLabelFormat("%.2f");
          axis_x->setTitleText("Water saturation");
          axis_x->setRange(0, 1);
          chart->addAxis(axis_x, Qt::AlignBottom);

          pc_axis_y_ = new QLogValueAxis;//static_cast<QLogValueAxis*>(chart_->axes(Qt::Vertical)[0]);
          pc_axis_y_->setLabelFormat("%.0e");
          pc_axis_y_->setTitleText("Capillary pressure, Pa");
          chart->addAxis(pc_axis_y_, Qt::AlignLeft);

          primary_.pc = new QLineSeries;
          primary_.pc->setName("primary");
          primary_.pc->clear();
          primary_.pc->setPointsVisible(true);

          secondary_.pc = new QLineSeries;
          secondary_.pc->setName("secondary");
          secondary_.pc->clear();
          secondary_.pc->setPointsVisible(true);

          chart->addSeries(primary_.pc);
          primary_.pc->attachAxis(axis_x);
          primary_.pc->attachAxis(pc_axis_y_);

          {
            auto* dummy = new QLineSeries;
            chart->addSeries(dummy);
            chart->addSeries(secondary_.pc);
            chart->removeSeries(dummy);
          }

          secondary_.pc->attachAxis(axis_x);
          secondary_.pc->attachAxis(pc_axis_y_);

          plot_tabs->addTab(pc_chart_view_, "Pc");
        }

        {
          kr_chart_view_ = new QChartView;
          kr_chart_view_->setRenderHint(QPainter::Antialiasing);
          kr_chart_view_->setBackgroundBrush(Qt::GlobalColor::white);

          {
            auto menu = std::make_unique<QMenu>();

            auto add = [&](auto caption, auto text_map) {
              connect(menu->addAction(caption), &QAction::triggered, [=] {
                QApplication::clipboard()->setText(QString::fromStdString(text_map()));
              });
            };

            static constexpr auto json_format =
              "[\n"
              "  {},\n"
              "  {}\n"
              "]\n";

            add("Copy (primary)", [this] { return model_.kr_to_plain(std::true_type{}); });
            add("Copy (secondary)", [this] { return model_.kr_to_plain(std::false_type{}); });

            menu->addSeparator();

            add("Copy JSON (primary)", [this] { return fmt::format(json_format, model_.kr_to_json(std::true_type{}), model_.kr_to_json(std::true_type{})); });
            add("Copy JSON (secondary)", [this] { return fmt::format(json_format, model_.kr_to_json(std::false_type{}), model_.kr_to_json(std::false_type{})); });
           
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

          primary_.kr0 = new QLineSeries;
          primary_.kr0->setPointsVisible(true);
          primary_.kr0->setName("primary-wat");

          primary_.kr1 = new QLineSeries;
          primary_.kr1->setPointsVisible(true);
          primary_.kr1->setName("primary-ph1");

          chart->addSeries(primary_.kr0);
          primary_.kr0->attachAxis(axis_x);
          primary_.kr0->attachAxis(kr_axis_y_);

          chart->addSeries(primary_.kr1);
          primary_.kr1->attachAxis(axis_x);
          primary_.kr1->attachAxis(kr_axis_y_);

          secondary_.kr0 = new QLineSeries;
          secondary_.kr0->setPointsVisible(true);
          secondary_.kr0->setName("secondary-wat");

          secondary_.kr1 = new QLineSeries;
          secondary_.kr1->setPointsVisible(true);
          secondary_.kr1->setName("secondary-ph1");

          chart->addSeries(secondary_.kr0);
          secondary_.kr0->attachAxis(axis_x);
          secondary_.kr0->attachAxis(kr_axis_y_);

          chart->addSeries(secondary_.kr1);
          secondary_.kr1->attachAxis(axis_x);
          secondary_.kr1->attachAxis(kr_axis_y_);

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


    



  public:
    void Init() {
      model_.init();
      
      InitGUI();

      InitPressure(model_.compute_pressure());

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

        for (const auto& p : model_.settings().primary.pc)
          darcy_pc_series->append(p.x(), p.y());
        darcy_pc_series->setName("microporous");
        darcy_pc_series->attachAxis(chart->axes(Qt::Horizontal)[0]);
        darcy_pc_series->attachAxis(pc_axis_y_);
        chart->legend()->markers(darcy_pc_series)[0]->setVisible(false);
        darcy_pc_series->setPen(QPen{Qt::gray, 1, Qt::DashLine});
      }

      using eq_tr = hydraulic_properties::equilateral_triangle;

      {
        using namespace std;

        auto min_r_cap_throat = numeric_limits<double>::max();

        for (size_t i{0}; i < model_.pni().pn().throat_count(); ++i)
          if (auto [l, r] = model_.pni().pn().throat_[attrib::adj][i]; model_.pni().pn().inner_node(r) && model_.pni().connected(l))
            min_r_cap_throat = min(min_r_cap_throat, model_.pni().pn().throat_[attrib::r_ins][i]);

        pc_axis_y_->setRange(
          // 1e4/*axis_y->min()*/,
          pow(10., floor(log10(1./eq_tr::r_cap_piston_with_films(0, ranges::max(model_.pni().pn().node_.span(attrib::r_ins)))))),
          pow(10., ceil(log10(max(
            1/(0.7*eq_tr::r_cap_piston_with_films(0, min_r_cap_throat)),
            model_.settings().primary.pc.empty() ? 0 : model_.settings().primary.pc.front().y()))))*1.01
        );
      }

      


      

      consumer_future_ = std::async(std::launch::async, [&] {
        model_.invasion_task().init();

        auto start = std::chrono::system_clock::now();

        auto pc_inv = model_.settings().primary.calc_pc_inv();

        auto invasion_future = std::async(std::launch::async, &invasion_task::launch_primary, &model_.invasion_task(),
          model_.absolute_rate(), model_.settings().theta, pc_inv);

        auto last_progress_idx = std::numeric_limits<idx1d_t>::max();

        auto update = [this, start, &pc_inv, &last_progress_idx] {
          // return;

          if (last_progress_idx == model_.invasion_task().progress_idx())
            return;

          last_progress_idx = model_.invasion_task().progress_idx();

          auto& state = model_.invasion_task().state();
          auto darcy_saturation = 1.0 - (pc_inv.empty() ? 0.0 : solve(pc_inv, 1/state.r_cap_global, dpl::extrapolant::flat));
          
          auto map_satur = [](double x) { return x/2. + 0.25; };



          // dpl::sfor<6>([&](auto face_idx) {
          //   dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
          //
          //   for (vtkIdType i = 0; auto idx1d : face.GetIndices()) {
          //     if (auto v_idx = voxel_t{idx1d};
          //       model_.pni().connected(v_idx) && state.config(model_.pni().net(v_idx)).phase() == phase_config::phase1())
          //       face.GetColorArray()->SetTypedComponent(i, 0, map_satur(darcy_saturation));
          //     else
          //       face.GetColorArray()->SetTypedComponent(i, 0, 0.0);
          //
          //     ++i;
          //   }
          // });
          //
          // using namespace attrib;
          //
          // for (macro_t i{0}; i < model_.pni().pn().node_count(); ++i)
          //   if (model_.pni().connected(i) && state.config(model_.pni().net(i)).phase() == phase_config::phase1())
          //     macro_colors->SetTypedComponent(*i, 0,
          //       map_satur(1.0 - eq_tr::area_corners(model_.settings().theta, state.r_cap(model_.pni().net(i)))/eq_tr::area(r_ins(model_.pni().pn(), i))));  // NOLINT(cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
          //   else
          //     macro_colors->SetTypedComponent(*i, 0, 0.0);
          //
          // {
          //   vtkIdType t_inner_idx = 0;
          //
          //   for (std::size_t i{0}; i < model_.pni().pn().throat_count(); ++i) {
          //     if (auto [l, r] = adj(model_.pni().pn(), i); model_.pni().pn().inner_node(r)) {
          //       if (model_.pni().connected(l) && state.config(i).phase() == phase_config::phase1())
          //         throat_colors->SetTypedComponent(t_inner_idx, 0,  // NOLINT(cppcoreguidelines-narrowing-conversions, clang-diagnostic-implicit-float-conversion)
          //           map_satur(1.0 - eq_tr::area_corners(model_.settings().theta, state.r_cap(i))/eq_tr::area(r_ins(model_.pni().pn(), i))));
          //       else
          //         throat_colors->SetTypedComponent(t_inner_idx, 0, 0.0);
          //
          //       ++t_inner_idx;
          //     }
          //   }
          // }






          
          
          QMetaObject::invokeMethod(this,
            [this, start] {
              using namespace std::chrono;

              UpdateStatus(model_.invasion_task().finished()
                ? fmt::format("done {}s", duration_cast<seconds>(system_clock::now() - start).count())
                : fmt::format("{:.1f} %", 100.*model_.invasion_task().progress_idx()/(*model_.pni().connected_count())));

              primary_.pc->clear();
              primary_.kr0->clear();
              primary_.kr1->clear();

              secondary_.pc->clear();
              secondary_.kr0->clear();
              secondary_.kr1->clear();

              for (auto p : model_.invasion_task().primary().pc)
                primary_.pc->append(p.x(), p.y());

              for (auto p : model_.invasion_task().secondary().pc)
                secondary_.pc->append(p.x(), p.y());

              for (auto [sw, kro, krg] : model_.invasion_task().primary().kr) {
                primary_.kr0->append(sw, kro);
                primary_.kr1->append(sw, krg);                
              }

              for (auto [sw, kro, krg] : model_.invasion_task().secondary().kr) {
                secondary_.kr0->append(sw, kro);
                secondary_.kr1->append(sw, krg);                
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

    void InitPressure(dpl::strong_vector<net_t, double> pressure) {
      // auto assembly = vtkSmartPointer<vtkAssembly>::New();

      {
        vtkSmartPointer<vtkActor> actor;

        std::tie(actor, macro_colors) = CreateNodeActor(pn(), lut_pressure_, 
          [&](macro_t i) {
            return pni().connected(i) ? /*0*/ pressure[pni().net(i)] : std::numeric_limits<double>::quiet_NaN();
          });

        macro_network_->AddPart(actor);
        
        std::tie(actor, throat_colors) = CreateThroatActor(pn(), lut_pressure_, [&](std::size_t i) {
          auto [l, r] = pn().throat_[attrib::adj][i];

          return pni().connected(l)
            ? /*0*/ (pressure[pni().net(l)] + pressure[pni().net(r)])/2
            : std::numeric_limits<double>::quiet_NaN();
        });

        macro_network_->AddPart(actor);
      }

      renderer_->AddActor(macro_network_);

      // std::ranges::fill(pressure, 0);

      {
        auto scale_factor = /*1.0*/pn().physical_size.x()/img().dim().x(); // needed for vtk 8.2 floating point arithmetics
            
        img_glyph_mapper_.Init(scale_factor);
        {
          dpl::strong_vector<voxel_t, bool> filter{img().size()};

          {
            idx3d_t ijk;
            auto& [i, j, k] = ijk;
            voxel_t idx1d{0};

            for (k = 0; k < img().dim().z(); ++k)
              for (j = 0; j < img().dim().y(); ++j)
                for (i = 0; i < img().dim().x(); ++i, ++idx1d) {
                  filter[idx1d] =
                    // pni().connected(voxel_t{idx1d}) &&
                    // !(i < img().dim().x()/2.*1.25 &&
                    //   j > img().dim().y()/2./1.25 &&
                    //   k > img().dim().z()/2./1.1) &&
                    img().dict.is_darcy(img().phase[idx1d]);      
                }
          }


          cout << "3D faces...";

          using seconds = std::chrono::seconds;
          using clock = std::chrono::high_resolution_clock;

          auto t0 = clock::now();

          img_glyph_mapper_.Populate(img().dim(), pn().physical_size/img().dim(),
            [&](idx1d_t idx1d) { return filter[voxel_t{idx1d}]; });


          std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());



          // img().dict


          // {

            auto* ptr = model_.settings().darcy.poro_perm.data();

          

            using poro_perm_t = runtime_settings::poro_perm_t;
            
            using namespace std::ranges;

            //
            // int k = 0;
            // for (auto pp : model_.settings().darcy.poro_perm) {
            //   if (!model_.settings().darcy.poro_perm[voxel_prop::phase_t(k)].is_nan()) 
            //     fmt::print("value: {}, poro: {}, perm: {}\n", k, pp.poro, pp.perm/presets::mD_to_m2);
            //   ++k;
            // }

            auto [minPERM, maxPERM] = minmax(
              subrange(ptr, ptr + 256) | views::filter([](poro_perm_t pp) { return !pp.is_nan(); }),
              {}, [](poro_perm_t x) { return x.perm; });

            // fmt::print("\n\nPERM MIN, MAX: ({} mD, {} mD)\n\n", minPERM.perm/presets::mD_to_m2, maxPERM.perm/presets::mD_to_m2);

            
          // }



          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);

            idx1d_t i = 0;
            for (auto idx1d : face.GetIndices())
              face.GetColorArray()->SetTypedComponent(i++, 0, 
                pni().connected(voxel_t{idx1d})
                  ?
                    /*0*/
                    pressure[pni().net(voxel_t{idx1d})]
                    
                    // (log10(model_.settings().darcy.perm(img().phase[voxel_t{idx1d}])) - log10(minPERM.perm))/(
                    //   
                    //   log10(maxPERM.perm) - log10(minPERM.perm))/1.25+0.1  /*/2+0.25*/

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
        0., pn().physical_size.x(),
        0., pn().physical_size.y(),
        0., pn().physical_size.z()};
    }
  };
}