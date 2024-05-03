#pragma once

#include "functions.h"
#include "modeller.hpp"

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
#include <QMenuBar>
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
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>
#include <vtkWindowToImageFilter.h>
#include <vtkProp3DCollection.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <algorithm>
#include <future>
#include <unordered_set>

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
    QSplitter* hsplit_;
    QVTKWidgetRef* qvtk_widget_;
    dpl::vtk::TidyAxes tidy_axes_;

    // vtkNew<vtkLookupTable> lut_velem_;
    // vtkNew<vtkLookupTable> lut_pore_solid_micro_;
    

    vtkNew<vtkAssembly> macro_network_;
    vtkFloatArray* macro_colors;
    vtkFloatArray* throat_colors;


    dpl::vtk::ImageDataGlyphMapper<idx1d_t> img_glyph_mapper_;





    modeller model_;

    pore_network_image& pni() { return model_.pni(); }
    pore_network& pn() { return pni().pn(); }
    image_data& img() { return pni().img(); }

    // pore_network pn_;
    // image_data img_;
    // pore_network_image pni_{pn_, img_};

    // runtime_settings settings_;


    std::array<double, 6> bounds_ = {0, 100, 0, 100, 0, 100};

    

    // std::string status_ = "<nothing>";
    // dpl::qt::property_editor::PropertyItem* status_property_item_;

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
    
    // void UpdateStatus(std::string text) {
    //   status_ = std::move(text);
    //   auto idx = tree_view_->model()->index(status_property_item_, 1);
    //   tree_view_->model()->dataChanged(idx, idx);
    // }

    std::future<void> consumer_future_;
    
    static void initLutNodeThroat(vtkLookupTable* lut) {
      lut->IndexedLookupOn();
      lut->SetNumberOfTableValues(2);
      lut->SetTableValue(0, 0.784314, 0.467419, 0.657556);
      lut->SetAnnotation(vtkVariant(0), "Node");
      lut->SetTableValue(1, 0.784314, 0.752624, 0.467419);
      lut->SetAnnotation(vtkVariant(1), "Throat");
    }

    static void initLutPhases(vtkLookupTable* lut, const parse::image_dict& d) {
      lut->IndexedLookupOn();  
      lut->SetNumberOfTableValues(256);
    
      for (int i = 0; i < 256; ++i) {
        lut->SetTableValue(i, 0.8, 0.8, 0.8);
        lut->SetAnnotation(vtkVariant(i), "Darcy");
      }
    
      // lut->SetTableValue(d.pore, 0.784314, 0.467419, 0.657556);
      
      lut->SetTableValue(d.pore, 194/255., 105/255., 158/255.);
      lut->SetAnnotation(vtkVariant(d.pore), "Void");
      lut->SetTableValue(d.solid, 0.5, 0.5, 0.5);
      lut->SetAnnotation(vtkVariant(d.solid), "Solid");
    }

    static void initLutPhases(vtkLookupTable* lut, const parse::image_dict& d,
      const dpl::strong_array<voxel_prop::phase_t, poro_perm_t>& poro_perms) {

      // lut->IndexedLookupOn();  
      lut->SetNumberOfTableValues(256);

      for (int i = 0; i < 256; ++i) {
        if (auto& pp = poro_perms[voxel_prop::phase_t(i)]; pp.is_nan())  // NOLINT(clang-diagnostic-implicit-int-conversion)
          lut->SetTableValue(i, 0.8, 0.8, 0.8);
        else
          lut->SetTableValue(i, pp.color.redF(), pp.color.greenF(), pp.color.blueF());

        // lut->SetAnnotation(vtkVariant(i), "Darcy");
      }

      lut->SetTableValue(d.pore, 0.784314, 0.467419, 0.657556);
      // lut->SetAnnotation(vtkVariant(d.pore), "Void");
      lut->SetTableValue(d.solid, 0.5, 0.5, 0.5);
      // lut->SetAnnotation(vtkVariant(d.solid), "Solid");

      lut->SetTableRange(0, 255);
    }

    



    static void InitLutVelems(vtkLookupTable* lut, std::int32_t max) {
      auto count = max - (-2)/**min*/ + 1;

      // lut->IndexedLookupOn();
      lut->SetNumberOfTableValues(count);
      
      for (int32_t i = 0; i < count; ++i) {
        auto coef = static_cast<double>(i*45%count)/count;
        auto color = QColor::fromHsl(coef*255, 175, 122);                     // NOLINT(cppcoreguidelines-narrowing-conversions)
        lut->SetTableValue(i, color.redF(), color.greenF(), color.blueF());   // NOLINT(clang-diagnostic-double-promotion)
      }
    
      lut->SetTableValue(0, 0.784314, 0.752624, 0.467419);      // -2 velem value, throat
      lut->SetTableValue(1, 0.8, 0.8, 0.8);                     // -1 velem value
      lut->SetTableRange(-2/**min*/, max);
    }


    class DisplayItem// : BaseWidgetItem<DisplayPropertyItem, std::string>
    {
      static inline auto fields_ = std::to_array<std::string>({ "one", "two", "three" });

    public:
      using ComboItem = std::string;

      Widget* widget;
      std::uint8_t index = 0;

      static constexpr auto Name() { return "Display"; }

      ComboItem* Get() const {
        return &fields_[index];
      }

      void Set(std::uint8_t i) {
        index = i;
        widget->render_window_->Render();
      }

      // void Set(ComboItem* prop) {
      //   std::cout << "\n\nPTR\n\n";
      //
      //   index = &*std::ranges::find(fields_, *prop) - fields_.data();
      //
      //   // widget->ScalarBar()->DrawColorBarOn();
      //   // widget->ScalarBar()->DrawAnnotationsOn();
      //   // widget->ScalarBar()->Modified();
      //
      //   // prop->SetActive();
      //   
      //   // widget->HardRefreshPropertyTree();
      //   widget->render_window_->Render();
      // }

      static auto Items() {
        std::vector<std::string*> items;

        for (auto& s : fields_)
          items.push_back(&s);

        return items;
      }

      static std::string Caption(ComboItem* item) {
        return *item;
      }
    };

    void initMainMenu() {
      {
        auto* cat = menuBar()->addMenu("File");

        connect(cat->addAction("Exit"), &QAction::triggered, [this] { close(); });
      }

      {
        auto* cat = menuBar()->addMenu("Edit");

        using namespace dpl::vtk;

        {
          auto act = addAction("Copy Camera");
          cat->addAction(act);
          act->setShortcut(QKeySequence{Qt::CTRL | Qt::SHIFT | Qt::Key_C});
          connect(act, &QAction::triggered, [this] {
            QApplication::clipboard()->setText(QString::fromStdString(
              nlohmann::json(CameraSettings{renderer_->GetActiveCamera()}).dump(2)));
          });
        }

        {
          auto act = addAction("Paste Camera");
          cat->addAction(act);
          act->setShortcut(QKeySequence{Qt::CTRL | Qt::SHIFT | Qt::Key_V});
          connect(act, &QAction::triggered, [this] {
            try {
              CameraSettings{
                nlohmann::json::parse(QApplication::clipboard()->text().toStdString())
              } >> renderer_->GetActiveCamera();
              renderer_->ResetCameraClippingRange();
              render_window_->Render();
            }
            catch (...) {}  // NOLINT(bugprone-empty-catch)
          });
        }

        cat->addSeparator();

        {
          auto* act = addAction("Save Screenshot");
          cat->addAction(act);
          act->setShortcut(QKeySequence{Qt::CTRL | Qt::SHIFT | Qt::Key_P});
          connect(act, &QAction::triggered, [this] {
            saveScreenshot();
          });
        }
      }

      {
        auto* cat = menuBar()->addMenu("Help");
      }
    }

    void saveScreenshot() {
      using namespace std::filesystem;

      create_directories("screenshots");

      int number = 0;

      if (number == 0) {
        std::regex reg(R"(screenshot_(\d+)\.png)", std::regex_constants::icase);
        std::smatch match;

        
        for (const auto& entry : directory_iterator("screenshots")) {
          auto stem = entry.path().filename().string();

          if (std::regex_search(stem, match, reg))
            number = std::max(number, std::stoi(match.str(1)));
        }
      }

      vtkNew<vtkWindowToImageFilter> filter;

      filter->SetInput(render_window_);
      filter->SetInputBufferTypeToRGB(); //also record the alpha (transparency) channel
      filter->ReadFrontBufferOff(); // read from the back buffer
      filter->Update();

      vtkNew<vtkPNGWriter> writer;

      
      writer->SetFileName(
        fmt::format("screenshots\\screenshot_{:03}.png", number + 1).c_str()
        // (boost::format("screenshots\\screenshot_%03i.png") % ++number).str().c_str()
      );
      writer->SetInputConnection(filter->GetOutputPort());
      writer->Write();
    }

    auto* viewport() const { return qvtk_widget_; }

    void setViewportWidth(int width) {
      auto delta = width - viewport()->width();
      
      auto sizes = hsplit_->sizes();
      sizes[0] -= delta;
      sizes[1] += delta;
      
      hsplit_->setSizes(sizes);
    }


    void keyPressEvent(QKeyEvent *event) override {
      if (event->key() == Qt::Key_F11) {
        if (isFullScreen()) {
          menuBar()->show();
          showNormal();
        }
        else {
          showFullScreen();
          menuBar()->hide();
        }
      }
      // else if (event->keyCombination() == QKeyCombination(Qt::CTRL | Qt::SHIFT, Qt::Key_P)           /*modifiers() == (Qt::ControlModifier | Qt::ShiftModifier) && event->key() == Qt::Key_P*/) {
      //   saveScreenshot();
      //   std::cout << "keyPressEvent";
      // }
      else
        QWidget::keyPressEvent(event);
    }

    void initGUI() {
      initMainMenu();

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

        using namespace dpl::qt::property_editor;

        model->emplace_back<ItemFunctor<bool>>(
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
          });

        model->emplace_back<ItemFunctor<bool>>(
          "Network",
          [this] {
            return static_cast<bool>(macro_network_->GetVisibility());
          },
          [this](bool v) {
            macro_network_->SetVisibility(v);
            render_window_->Render();
          });

        model->emplace_back<ItemFunctor<bool>>(
          "Darcy",
          [this] {
            return static_cast<bool>(std::get<0>(img_glyph_mapper_.faces_).GetActor()->GetVisibility());
          },
          [this](bool v) {
            dpl::sfor<6>([this, v](auto i) {
              static_cast<dpl::vtk::GlyphMapperFace<idx1d_t>&>(std::get<i>(img_glyph_mapper_.faces_))
                .GetActor()->SetVisibility(v);
              render_window_->Render();
            });
          });

        model->emplace_back<ItemFunctor<bool>>(
          "Axes",
          [this] {
            return tidy_axes_.GetVisibility();
          },
          [this](bool v) {
            tidy_axes_.SetVisibility(v);
            render_window_->Render();
          });

        model->emplace_back<ItemFunctor<double>>(
          "Camera angle",
          [this] {
            return renderer_->GetActiveCamera()->GetViewAngle();
          },
          [this](double v) {
            renderer_->GetActiveCamera()->SetViewAngle(v);
            render_window_->Render();
          });

        

        // model->emplace_back<ItemFunctor<bool>>(
        //   "Main menu",
        //   [this] {
        //     return menuBar()->isVisible();
        //   },
        //   [this](bool v) {
        //     if (v)
        //       menuBar()->show();
        //     else
        //       menuBar()->hide();
        //   });

        model->emplace_back<ItemFunctor<dpl::vector2i>>(
          "Viewport, px",
          [this] { return dpl::vector2i{viewport()->width(), viewport()->height()}; },
          [this](const dpl::vector2i& v) {
            if (isMaximized() || isFullScreen())
              setViewportWidth(v.x());
            else {
              QSize desired{
                width() + v.x() - viewport()->width(),
                height() + v.y() - viewport()->height()
              };

              setGeometry(
                QStyle::alignedRect(
                  Qt::LeftToRight,
                  Qt::AlignLeft | Qt::AlignTop,
                  desired,
                  geometry()));
            }
          });




  //       struct ViewportItem : BaseWidgetItem<ViewportItem, dpl::vector2i>
  // {
  //   constexpr auto Name() { return "Viewport, px"; }
  //
  //   Type Get_() { return {widget->ViewportWidth(), widget->ViewportHeight()}; }
  //   void Set_(const Type& v) {
  //     if (widget->isMaximized())
  //       widget->SetOpenGLWidgetWidth(v.x());
  //     else {
  //       QSize desired = {
  //         widget->width() + v.x() - widget->ViewportWidth(),
  //         widget->height() + v.y() - widget->ViewportHeight()
  //       };
  //
  //       widget->setGeometry(
  //         QStyle::alignedRect(
  //           Qt::LeftToRight,
  //           Qt::AlignLeft | Qt::AlignTop, desired,
  //           widget->geometry()));
  //     }
  //     
  //   }
  // };



        // model->AddItem(DisplayItem{this});







        // int k =  3;
        // di.Set(k);
        // Settable<DisplayItem, std::size_t>;
        // constexpr auto q = Settable<DisplayItem, int>;


        // constexpr auto foo = ;



        // using foo = decltype(static_cast<void (DisplayItem::*)(std::string)>(&DisplayItem::Set));

        // static_cast<void (DisplayItem::*)(int)>(&DisplayItem::Get)

        // auto q = decltype(&DisplayItem::Get);

        // constexpr auto qq = std::is_invocable_v<foo, const DisplayItem&, std::string>;





        // status_property_item_ = model->AddItem(
        //   dpl::qt::property_editor::ItemFunctor<std::string>{"Status", [this] { return status_; }});




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

            add("Copy JSON (primary)", [this] { return fmt::format(
              json_format,
              model_.kr_to_json(std::true_type{}, std::false_type{}),
              model_.kr_to_json(std::true_type{}, std::true_type{})); });

            add("Copy JSON (secondary)", [this] { return fmt::format(
              json_format,
              model_.kr_to_json(std::false_type{}, std::false_type{}),
              model_.kr_to_json(std::false_type{}, std::true_type{})); });
           
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

        splitter hsplit{dir::horizontal};
        hsplit_ = hsplit;
        
        setCentralWidget(
          hsplit
          << (
            splitter{dir::vertical}
            << tree_view_ << stretch{1}
            << plot_tabs/*chart_view_*/ << stretch{0}) << stretch{0}
          << qvtk_widget_ << stretch{1}
        );
      }
    }


    
  public:
    void Init(const std::filesystem::path& input) {
      model_.init(input);
      
      initGUI();

      

      // dpl::vtk::PopulateLutRedWhiteBlue(generic_lut_);
      dpl::strong_vector<net_t, double> pressure;

      model_.prepare();

      if (const auto& name = model_.settings().report.display; name == "velem") {
        vtkNew<vtkLookupTable> lut;

        InitLutVelems(lut, *pn().node_count() + 1);  // NOLINT(cppcoreguidelines-narrowing-conversions)

        velem_map map{&img()};

        InitNetworkImage(
          [this](voxel_t v) { return img().dict.is_darcy(img().phase[v]); },
          map,
          lut.Get(),
          lut.Get());

        macro_network_->RemovePart(macro_network_->GetParts()->GetLastProp3D());
      }
      else if (name == "phases") {
        vtkNew<vtkLookupTable> lut_image, lut_network;

        initLutPhases(lut_image, img().dict, model_.settings().darcy.poro_perm);
        initLutNodeThroat(lut_network);

        phases_map map{&img()};

        InitNetworkImage(
          [this](voxel_t v) { return img().dict.is_darcy(img().phase[v]); },
          map,
          lut_image,
          lut_network.Get());
      }
      else if (name == "phases-triplet") {
        vtkNew<vtkLookupTable> lut_image, lut_network;

        initLutPhases(lut_image, img().dict);
        initLutNodeThroat(lut_network);

        phases_map map{&img()};

        InitNetworkImage(
          [this](voxel_t) { return true; },
          // [this](voxel_t v) { return img().dict.is_darcy(img().phase[v]); },
          map,
          lut_image.Get(),
          lut_network.Get());
      }
      else if (name == "permeability") {
        vtkNew<vtkLookupTable> lut;
        dpl::vtk::PopulateLutRedWhiteBlue(lut);
        lut->SetTableRange(0, 1);

        permeability_map map{img(), model_.settings().darcy.poro_perm};

        InitNetworkImage(
          [this](voxel_t v) { return img().dict.is_darcy(img().phase[v]); },
          map,
          lut.Get(),
          lut.Get());
        
        //
        // int k = 0;
        // for (auto pp : model_.settings().darcy.poro_perm) {
        //   if (!model_.settings().darcy.poro_perm[voxel_prop::phase_t(k)].is_nan()) 
        //     fmt::print("value: {}, poro: {}, perm: {}\n", k, pp.poro, pp.perm/presets::mD_to_m2);
        //   ++k;
        // }

        // std::span(ptr, 256)// TODO

        // fmt::print("\n\nPERM MIN, MAX: ({} mD, {} mD)\n\n", minPERM.perm/presets::mD_to_m2, maxPERM.perm/presets::mD_to_m2);
      }
      else {
        vtkNew<vtkLookupTable> lut;
        dpl::vtk::PopulateLutRedWhiteBlue(lut);
        lut->SetTableRange(0, 1);

        pressure = model_.compute_pressure();

        pressure_map map{&pni(), &pressure};

        InitNetworkImage(
          [this](voxel_t v) { return img().dict.is_darcy(img().phase[v]); },
          map,
          lut.Get(),
          lut.Get());
      }

      

      renderer_->ResetCamera(bounds_.data());
      renderer_->GetActiveCamera()->Zoom(0.70);
      // renderer_->ResetCameraClippingRange();
      
      tidy_axes_.Init(renderer_.Get());

      connect(qvtk_widget_, &QVTKWidgetRef::resized, this, [this] { tidy_axes_.RefreshAxes(); });

      // tidy_axes_.SetScale(1.e-6/*startup.image.resolution*/);
      // tidy_axes_.SetFormat(".2e");
      tidy_axes_.Build(bounds_.data());

      // fmt::print("CAMERA ANGLE: {}", renderer_->GetActiveCamera()->GetViewAngle());


      // {
      //   net_idx_t connected = pni_.connected_total_count();
      //
      //   occupancy_arrays_.macro.resize(*connected);
      //   for (net_idx_t i{0}; i < connected; ++i)
      //     occupancy_arrays_.macro[i] = 0;
      // }

      if (model_.settings().report.invasion_percolation) {
        if (!pressure)
          model_.compute_pressure();

        LaunchInvasion();
      }
      else {
        pc_axis_y_->setRange(1, 1000);
      }
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
            model_.settings().primary.pc ? model_.settings().primary.pc.front().y() : 0))))*1.01
        );
      }

      


      

      consumer_future_ = std::async(std::launch::async, [&] {
        model_.get_invasion_task().init();

        auto start = std::chrono::system_clock::now();

        // auto pc_inv = model_.settings().primary.calc_pc_inv();

        auto invasion_future = std::async(std::launch::async, &invasion_task::launch_primary, &model_.get_invasion_task(),
          model_.absolute_rate(),
          model_.settings().theta,
          model_.settings().primary.pc.inverse_unique());

        auto last_progress_idx = std::numeric_limits<idx1d_t>::max();

        auto update = [this, start, &last_progress_idx] {
          // return;

          if (last_progress_idx == model_.get_invasion_task().progress_idx())
            return;

          last_progress_idx = model_.get_invasion_task().progress_idx();

          
          
          

          

          if (model_.settings().report.display == "saturation") {
            const auto& state = model_.get_invasion_task().state();
            const auto& pni = model_.pni();
            
            static constexpr auto map_satur = [](double x) -> float { return x/2 + 0.25; };  // NOLINT(clang-diagnostic-implicit-float-conversion)

            auto pc_inv = (
              model_.get_invasion_task().primary_finished()
                ? model_.settings().secondary
                : model_.settings().primary).pc.inverse_unique();

            dpl::sfor<6>([&](auto face_idx) {
              dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
            
              for (vtkIdType i = 0; auto idx1d : face.GetIndices()) {
                auto sw = 0.0;

                if (voxel_t voxel{idx1d}; pni.connected(voxel))
                  if (auto net = pni.net(voxel); state.config(net).phase() == phase_config::phase1())
                    sw = 1.0 - pc_inv.solve(1/state.r_cap(net));

                face.GetColorArray()->SetTypedComponent(i++, 0, map_satur(sw));
              }
            });
            
            using namespace attrib;

            auto& pn = pni.pn();

            for (macro_t m{0}; m < pn.node_count(); ++m) {
              auto sw = 0.0;

              if (pni.connected(m))
                if (auto net = pni.net(m); state.config(net).phase() == phase_config::phase1())
                  sw = 1.0 - eq_tr::area_corners(model_.settings().theta, state.r_cap(pni.net(m)))/eq_tr::area(r_ins(pn, m));

              macro_colors->SetTypedComponent(*m, 0, map_satur(sw));
            }

            for (vtkIdType i = 0; throat_t t : pn.throats()) // for (throat_t t{0}; t < pn.throat_count(); ++t)
              if (auto [l, r] = adj(pn, t); pn.inner_node(r)) {
                auto sw = 0.0;

                if (pni.connected(l) && state.config(t).phase() == phase_config::phase1())
                  sw = 1.0 - eq_tr::area_corners(model_.settings().theta, state.r_cap(t))/eq_tr::area(r_ins(pn, t));

                throat_colors->SetTypedComponent(i++, 0, map_satur(sw));
              }
          }





          
          
          QMetaObject::invokeMethod(this,
            [this, start] {
              using namespace std::chrono;

              // UpdateStatus(model_.invasion_task().finished()
              //   ? fmt::format("done {}s", duration_cast<seconds>(system_clock::now() - start).count())
              //   : fmt::format("{:.1f} %", 100.*model_.invasion_task().progress_idx()/(*model_.pni().connected_count())));

              primary_.pc->clear();
              primary_.kr0->clear();
              primary_.kr1->clear();

              secondary_.pc->clear();
              secondary_.kr0->clear();
              secondary_.kr1->clear();

              for (auto p : model_.get_invasion_task().primary().pc)
                primary_.pc->append(p.x(), p.y());

              for (auto p : model_.get_invasion_task().secondary().pc)
                secondary_.pc->append(p.x(), p.y());

              for (auto [sw, kro, krg] : model_.get_invasion_task().primary().kr) {
                primary_.kr0->append(sw, kro);
                primary_.kr1->append(sw, krg);                
              }

              for (auto [sw, kro, krg] : model_.get_invasion_task().secondary().kr) {
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


    struct velem_map
    {
      image_data* img;

      auto operator()(voxel_t i) const {
        auto velem = img->velem[i];
        return velem ? *velem : -1;
      }

      auto operator()(macro_t i) const {
        return *i;
      }

      auto operator()(throat_t) const {
        return -2;
      }
    };

    struct phases_map
    {
      image_data* img;

      auto operator()(voxel_t i) const {
        return *img->phase[i];
      }

      auto operator()(macro_t) const {
        return 0;
      }

      auto operator()(throat_t) const {
        return 1;
      }
    };

    struct pressure_map
    {
      pore_network_image* pni;
      dpl::strong_vector<net_t, HYPRE_Real>* data;

      auto operator()(macro_voxel_t auto i) const {
        return pni->connected(i) ? (*data)[pni->net(i)] : std::numeric_limits<HYPRE_Real>::quiet_NaN();
      }

      auto operator()(throat_t i) const {
        auto [l, r] = attrib::adj(pni->pn(), i);

        return pni->connected(l)
          ? ((*data)[pni->net(l)] + (*data)[pni->net(r)])/2
          : std::numeric_limits<HYPRE_Real>::quiet_NaN();
      }

      // auto operator()(voxel_t i) const {
      //   return pni->connected(i) ? (*pressure)[pni->net(i)] : std::numeric_limits<double>::quiet_NaN();
      //
      //   pni().connected(voxel_t{idx1d})
      //         ?
      //           /*0*/
      //           pressure[pni().net(voxel_t{idx1d})]
      //           
      //           // (log10(model_.settings().darcy.perm(img().phase[voxel_t{idx1d}])) - log10(minPERM.perm))/(
      //           //   
      //           //   log10(maxPERM.perm) - log10(minPERM.perm))/1.25+0.1  /*/2+0.25*/
      //
      //         : std::numeric_limits<double>::quiet_NaN()
      // }
    };


    struct permeability_map
    {
    private:
      const image_data* img;
      const dpl::strong_array<voxel_prop::phase_t, poro_perm_t>* poro_perm;
      double log10_min;
      double log10_diff; 
    
    public:
      permeability_map(
        const image_data& img,
        const dpl::strong_array<voxel_prop::phase_t, poro_perm_t>& poro_perm)
        : img(&img),
          poro_perm(&poro_perm) {
    
        auto* ptr = poro_perm.data();
            
        using namespace std::ranges;
    
        if (auto range = subrange(ptr, ptr + 256) | views::filter([](const poro_perm_t& pp) { return !pp.is_nan(); }); !empty(range)) {
          auto [min, max] = minmax(range, {}, [](const poro_perm_t& x) { return x.perm; });
          log10_min = log10(min.perm);
          log10_diff = log10(max.perm) - log10_min;
        }
        // else {
        //   log10_min = std::numeric_limits<double>::quiet_NaN();
        //   log10_diff = std::numeric_limits<double>::quiet_NaN();
        // }
      }
    
    
      auto operator()(macro_throat_t auto) const {
        return std::numeric_limits<double>::quiet_NaN();
      }
    
      auto operator()(voxel_t i) const {
        return (log10((*poro_perm)[img->phase[i]].perm) - log10_min)/log10_diff/1.25+0.1;  /*/2+0.25*/
      }
    };


    void InitNetworkImage(const auto& filter, const auto& map, vtkLookupTable* lut_img, vtkLookupTable* lut_network) { /* dpl::strong_vector<net_t, double> pressure*/
      {
        vtkSmartPointer<vtkActor> actor;

        std::tie(actor, macro_colors) = xpm::CreateNodeActor(pn(), lut_network, map);
        macro_network_->AddPart(actor);

        std::tie(actor, throat_colors) = xpm::CreateThroatActor(pn(), lut_network, map);
        macro_network_->AddPart(actor);
      }

      renderer_->AddActor(macro_network_);

      {
        img_glyph_mapper_.Init(
          /*
           * scale_factor
           *
           * needed for vtk 8.2 floating point arithmetics
           */
          /*1.0*/pn().physical_size.x()/img().dim().x()  // NOLINT(clang-diagnostic-implicit-int-float-conversion)
        );


        {
          dpl::strong_vector<voxel_t, bool> filter_cache{img().size()};

          {
            idx3d_t ijk;
            auto& [i, j, k] = ijk;
            voxel_t idx1d{0};

            for (k = 0; k < img().dim().z(); ++k)
              for (j = 0; j < img().dim().y(); ++j)
                for (i = 0; i < img().dim().x(); ++i, ++idx1d) {
                  filter_cache[idx1d] = filter(idx1d);
                    // pni().connected(voxel_t{idx1d}) &&
                    // !(i < img().dim().x()/2.*1.25 &&
                    //   j > img().dim().y()/2./1.25 &&
                    //   k > img().dim().z()/2./1.1) &&
                }
          }


          cout << "3D faces...";

          using seconds = std::chrono::seconds;
          using clock = std::chrono::high_resolution_clock;

          auto t0 = clock::now();

          img_glyph_mapper_.Populate(img().dim(), pn().physical_size/img().dim(), [&](idx1d_t i) { return filter_cache[voxel_t{i}]; });

          std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());

          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);

            for (vtkIdType i = 0; auto idx1d : face.GetIndices())
              face.GetColorArray()->SetTypedComponent(i++, 0, map(voxel_t{idx1d}));

            auto* glyphs = face.GetGlyphMapper();
            auto* actor = face.GetActor();

            glyphs->SetLookupTable(lut_img);
              
            glyphs->SetColorModeToMapScalars();
            glyphs->UseLookupTableScalarRangeOn();
              
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
      
      tidy_axes_.SetScale(1.e-6);

      bounds_ = {
        0., pn().physical_size.x(),
        0., pn().physical_size.y(),
        0., pn().physical_size.z()};
    }
  };
}