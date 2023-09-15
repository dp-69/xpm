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



#undef LoadImage



namespace xpm
{
  


  

  


  
  #if (VTK_MAJOR_VERSION == 8)
    using QVTKWidgetRef = QVTKOpenGLWidget;
  #elif (VTK_MAJOR_VERSION == 9)
    using QVTKWidgetRef = QVTKOpenGLNativeWidget;
  #endif


  //   template<typename Index = std::size_t>
  // class graph_generator
  // {
  //   dc_graph g_;
  //   std::unique_ptr<std::size_t[]> shift_;
  //   std::size_t edge_count_ = 0;
  //
  //   using traits = strong_traits<Index>;
  //
  // public:
  //   explicit graph_generator(std::size_t vertex_count)
  //     : g_{vertex_count} {
  //
  //     shift_= std::make_unique<std::size_t[]>(vertex_count + 1);
  //     std::fill_n(shift_.get(), vertex_count + 1, 0);
  //   }
  //
  //   void reserve(Index l, Index r) {
  //     ++shift_[traits::get(l) + 1]; 
  //     ++shift_[traits::get(r) + 1];
  //     edge_count_ += 1;
  //   }
  //
  //   void allocate() {
  //     for (auto i = 0; i < g_.vertex_count(); ++i) {
  //       shift_[i + 1] += shift_[i];
  //       g_.directed_edges_end()[i] = shift_[i + 1];
  //     }
  //
  //     g_.allocate_edges(edge_count_);
  //   }
  //
  //   void set(Index l, Index r) {
  //     auto set_pair = [](vertex* v, vertex* u, directed_edge* vu, directed_edge* uv) {            
  //       vu->v1 = u;
  //       uv->v1 = v;
  //
  //       vu->opposite = uv;
  //       uv->opposite = vu;
  //     };
  //
  //     set_pair(
  //       g_.get_vertex(traits::get(l)),
  //       g_.get_vertex(traits::get(r)),
  //       g_.get_directed_edge(shift_[traits::get(l)]++),
  //       g_.get_directed_edge(shift_[traits::get(r)]++));
  //   }
  //
  //   auto acquire() {
  //     return std::move(g_);
  //   }
  // };

  // class macro_generator
  // {
  //   idx1d_t macro_count_;
  //   std::unique_ptr<std::size_t[]> shift_;
  //   std::size_t edge_count_ = 0;
  //
  //   using traits = dpl::strong_traits<macro_idx>;
  //
  // public:
  //   explicit macro_generator(idx1d_t macro_count)
  //     : macro_count_{macro_count}/*: g_{vertex_count}*/ {
  //
  //     shift_= std::make_unique<std::size_t[]>(macro_count + 1);
  //     std::fill_n(shift_.get(), macro_count + 1, 0);
  //   }
  //
  //   void reserve(macro_idx l, macro_idx r) {
  //     ++shift_[*l + 1]; 
  //     ++shift_[*r + 1];
  //     edge_count_ += 1;
  //   }
  //
  //   auto allocate_and_acquire() {
  //     for (auto i = 0; i < macro_count_; ++i)
  //       shift_[i + 1] += shift_[i];
  //
  //     return std::move(shift_);
  //   }
  // };

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
    vtkNew<vtkLookupTable> lut_velem_;
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

    startup_settings startup_;

    invasion_task invasion_task_{pni_};


    std::string status_ = "<nothing>";
    dpl::qt::property_editor::PropertyItem* status_property_item_;

    QLineSeries* sweep_series_;
    
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

      renderer_->SetBackground(v3d{1});


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




        auto sweep_chart_view_ = new QChartView;
        sweep_chart_view_->setRenderHint(QPainter::Antialiasing);
        sweep_chart_view_->setBackgroundBrush(Qt::GlobalColor::white);

        sweep_series_ = new QLineSeries;

        auto sweep_chart_ = new QChart;
        // sweep_chart_->legend()->hide();
        sweep_chart_->addSeries(sweep_series_);
        sweep_chart_->createDefaultAxes();

        sweep_chart_view_->setChart(sweep_chart_);

        auto* axis_x = static_cast<QValueAxis*>(sweep_chart_->axes(Qt::Horizontal)[0]);
        // axis_x->setLabelsFont(scaled_font);
        axis_x->setLabelFormat("%.2f");
        axis_x->setTitleText("Water saturation");
        // axis_x->setTitleFont(font);
        axis_x->setRange(0, 1);
        // axis_x->setLabelsBrush(black_brush);
        // axis_x->setTitleBrush(black_brush);


        auto* axis_y = static_cast<QValueAxis*>(sweep_chart_->axes(Qt::Vertical)[0]);
        // axis_y->setLabelsFont(scaled_font);
        axis_y->setLabelFormat("%.0f");
        axis_y->setTitleText("Capillary pressure, Pa");
        // axis_y->setTitleFont(font);
        
        axis_y->setRange(0, 400000);
        // axis_y->setLabelsBrush(black_brush);
        // axis_y->setTitleBrush(black_brush);

        sweep_series_->clear();
        // sweep_series_->append(1, 0);

        // sweep_series_->append(0.95, 100);
        // sweep_series_->append(0.90, 500);
        // sweep_series_->append(0.85, 1000);
        // sweep_series_->append(0.80, 5000);

        sweep_series_->setPointsVisible(true);



        using namespace dpl::qt::layout;



        setCentralWidget(
          splitter{dir::horizontal}
          << (
            splitter{dir::vertical}
            << tree_view_ << stretch{1}
            << sweep_chart_view_ << stretch{0}) << stretch{0}
          << qvtk_widget_ << stretch{1}
        );
      }
    }



    auto ProcessImage(const std::filesystem::path& p) {
      namespace fs = std::filesystem;

      auto filename = p.filename();

      if (auto copy_fn = "pnextract"/filename; fs::absolute(copy_fn) != fs::absolute(p))
        copy(p, copy_fn, fs::copy_options::update_existing);

      auto files = {"_link1.dat", "_link2.dat", "_node1.dat", "_node2.dat", "_VElems.raw"};

      auto prev = fs::current_path();
      fs::current_path(prev/"pnextract");
      auto network_dir = p.stem();

      if (std::ranges::all_of(files, [&](const auto& file) { return fs::exists(network_dir/file); })) {
        std::cout << "using cached network\n\n";
      }
      else {
        std::cout << "=========== pnextract's network extraction begin ===========\n";

        /*auto value = */std::system( // NOLINT(concurrency-mt-unsafe)
          fmt::format("{} {}", fs::current_path()/"pnextract", filename).c_str());
      
        fs::create_directory(network_dir);
        for (fs::path f : files)
          fs::rename(f, network_dir/f);

        fs::remove("_VElems.mhd");

        std::cout << "=========== pnextract's network extraction end ===========\n\n";
      }

      network_dir = fs::absolute(network_dir);

      fs::current_path(prev);

      return network_dir;
    }


  public:
    void Init() {
      std::locale::global(std::locale("en_US.UTF-8"));

      if (std::filesystem::exists("config.json"))
        startup_.load(nlohmann::json::parse(std::ifstream{"config.json"}, nullptr, true, true));

      InitGUI();

      ComputePressure();

      renderer_->ResetCamera(bounds_.data());
      renderer_->GetActiveCamera()->Zoom(0.70);
      // renderer_->ResetCameraClippingRange();
      
      tidy_axes_.Init(renderer_.Get());

      connect(qvtk_widget_, &QVTKWidgetRef::resized, this, [this]() { tidy_axes_.RefreshAxes(); });

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
      invasion_task_.init();

      consumer_future_ = std::async(std::launch::async, [this] {
        auto start = std::chrono::system_clock::now();

        auto invasion_future = std::async(std::launch::async,
          &invasion_task::launch, &invasion_task_, startup_.microporous_poro);

        auto update = [this, start] {
          auto invaded = invasion_task_.invaded_array();

          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
            
            for (vtkIdType i = 0; auto idx1d : face.GetIndices()) {
              if (pni_.connected(voxel_idx_t{idx1d}) && invaded[pni_.net(voxel_idx_t{idx1d})])
                face.GetColorArray()->SetTypedComponent(i, 0, 0.5);
              ++i;
            }
          });
          
          for (macro_idx_t i{0}; i < pn_.node_count(); ++i)
            if (pni_.connected(i) && invaded[pni_.net(macro_idx_t{i})])
              macro_colors->SetTypedComponent(*i, 0,
                // 1 - props::area_of_films(theta, r_cap)/props::area(pn_.node_[attribs::r_ins][*i])
                0.5
              );
          
          for (vtkIdType throat_net_idx = 0; auto [l, r] : pn_.throat_.range(attribs::adj))
            if (pn_.inner_node(r)) {
              if (pni_.connected(l)
                && invaded[pni_.net(l)]
                && invaded[pni_.net(r)])
                throat_colors->SetTypedComponent(throat_net_idx, 0, 0.5);
          
              ++throat_net_idx;
            }
          
          QMetaObject::invokeMethod(this,
            [this, start] {
              using namespace std::chrono;

              UpdateStatus(invasion_task_.finished()
                ? fmt::format("done {}s", duration_cast<seconds>(system_clock::now() - start).count())
                : fmt::format("{:.1f} %", 100.*invasion_task_.inv_idx()/(*pni_.connected_count())));

              sweep_series_->clear();
              
              for (auto p : invasion_task_.pc_curve())
                sweep_series_->append(p.x(), p.y());

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




      // using clock = std::chrono::high_resolution_clock;
      // using seconds = std::chrono::seconds;
      //
      // using namespace dpl::graph;
      //
      // auto t0 = clock::now();
      // std::cout << "graph...";
      // dc_graph_ = pni_.generate_dc_graph<true>();
      // std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
      //   duration_cast<seconds>(clock::now() - t0).count(),
      //   dc_graph_.vertex_count(),
      //   dc_graph_.edge_count());
      //
      // auto t1 = clock::now();
      // std::cout << "euler tour...";
      // dc_context_.init_with_dfs(dc_graph_, dc_properties{dc_graph_});
      // std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t1).count());
      //
      //
      // std::cout << "full decremental connectivity... (async)\n";
      //
      // invasion_future_ = std::async(std::launch::async, [this] {
      //   auto t2 = clock::now();
      //
      //   auto theta = 0.0;
      //
      //   auto index_count = *pni_.connected_count();
      //
      //   using props = hydraulic_properties::equilateral_triangle_properties;
      //
      //   struct processed_t
      //   {
      //     dpl::strong_array<net_tag, bool> invaded_macro_voxel;
      //     dpl::strong_array<net_tag, bool> explored_macro_voxel;
      //
      //     processed_t(net_idx_t connected_count)
      //       : invaded_macro_voxel(connected_count), explored_macro_voxel(connected_count) {}
      //   };
      //
      //   double darcy_r_cap_const = 
      //     props::r_cap_piston_with_films(theta, std::ranges::min(pn_.node_.range(attribs::r_ins)))*0.95;
      //
      //   auto processed = std::make_shared<processed_t>(pni_.connected_count());
      //
      //   displ_queue queue;
      //
      //   for (std::size_t i{0}; i < pn_.throat_count(); ++i)
      //     if (auto [l, r] = pn_.throat_[attribs::adj][i]; pn_.inlet() == r) // inlet macro nodes TODO: voxels inlet needed
      //       if (pni_.connected(l)) {
      //         queue.insert(
      //           displ_elem::macro, *l, props::r_cap_piston_with_films(theta, pn_.node_[attribs::r_ins][*l]));
      //          
      //         processed->explored_macro_voxel[pni_.net(l)] = true;
      //       }
      //
      //   auto update_3d = [this, processed/*, theta*/](double r_cap) {
      //     dpl::sfor<6>([&](auto face_idx) {
      //       dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);
      //       
      //       for (vtkIdType i = 0; auto idx1d : face.GetIndices()) {
      //         if (pni_.connected(voxel_idx_t{idx1d}) && processed->invaded_macro_voxel[pni_.net(voxel_idx_t{idx1d})])
      //           face.GetColorArray()->SetTypedComponent(i, 0, 0.5);
      //         ++i;
      //       }
      //     });
      //
      //     for (macro_idx_t i{0}; i < pn_.node_count(); ++i)
      //       if (pni_.connected(i) && processed->invaded_macro_voxel[pni_.net(macro_idx_t{i})])
      //         macro_colors->SetTypedComponent(*i, 0,
      //           // 1 - props::area_of_films(theta, r_cap)/props::area(pn_.node_[attribs::r_ins][*i])
      //           0.5
      //         );
      //
      //     for (vtkIdType throat_net_idx = 0; auto [l, r] : pn_.throat_.range(attribs::adj))
      //       if (pn_.inner_node(r)) {
      //         if (pni_.connected(l)
      //           && processed->invaded_macro_voxel[pni_.net(l)]
      //           && processed->invaded_macro_voxel[pni_.net(r)])
      //           throat_colors->SetTypedComponent(throat_net_idx, 0, 0.5);
      //
      //         ++throat_net_idx;
      //       }
      //
      //     QMetaObject::invokeMethod(this,
      //       [this] {
      //         dpl::sfor<6>([&](auto face_idx) {
      //           std::get<face_idx>(img_glyph_mapper_.faces_).GetColorArray()->Modified();
      //         });
      //         throat_colors->Modified();
      //         macro_colors->Modified();
      //
      //         render_window_->Render();
      //       });
      //   };
      //
      //   std::future<void> update_future;
      //
      //   dc_properties dc_props{dc_graph_};
      //
      //   net_idx_t outlet_idx{index_count};
      //   et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);
      //
      //   constexpr auto delay = std::chrono::milliseconds{250};
      //   auto last = clock::now() - delay;
      //
      //   v3d inv_volume_coefs{0};
      //
      //   auto last_r_cap = queue.front().radius_cap;
      //
      //   auto total_pore_volume = 0.0;
      //
      //   for (macro_idx_t i{0}; i < pn_.node_count(); ++i)
      //     if (pni_.connected(i))
      //       total_pore_volume += pn_.node_[attribs::volume][*i];
      //
      //   auto unit_darcy_pore_volume = (pn_.physical_size/img_.dim).prod()*startup_.microporous_poro;
      //
      //   for (voxel_idx_t i{0}; i < img_.size; ++i)
      //     if (pni_.connected(i))
      //       total_pore_volume += unit_darcy_pore_volume;
      //
      //   dpl::vector2d last_pc_point{2, -1};
      //
      //   bool darcy_invaded = false;
      //
      //   auto add_to_plot = [this](dpl::vector2d p) {
      //     QMetaObject::invokeMethod(this, [=, this] {
      //       if (const auto& pts = sweep_series_->points(); pts.size() > 1
      //         && std::abs(pts[pts.size() - 1].y() - p.y()) < 1e-6
      //         && std::abs(pts[pts.size() - 2].y() - p.y()) < 1e-6)
      //         sweep_series_->replace(pts.size() - 1, p.x(), p.y());  // NOLINT(cppcoreguidelines-narrowing-conversions)
      //       else
      //         sweep_series_->append(p.x(), p.y());
      //     });
      //   };
      //
      //   for (idx1d_t inv_idx = 0; !queue.empty(); ++inv_idx) {
      //     auto inv_volume = inv_volume_coefs[0] + inv_volume_coefs[2]*last_r_cap*last_r_cap;
      //
      //     if (dpl::vector2d pc_point{1 - inv_volume/total_pore_volume, 1/last_r_cap};
      //       std::abs(last_pc_point.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point.y() - pc_point.y()) > 1e4) {
      //       last_pc_point = pc_point;
      //       add_to_plot(pc_point);
      //     }
      //
      //
      //     if (darcy_invaded) {
      //       if (inv_idx % ((index_count - 1)/50) == 0) {
      //         // if (reference_count == pni_.connected_macro_count())
      //         //   if (auto diff = duration_cast<seconds>(clock::now() - last); diff < delay)
      //         //     std::this_thread::sleep_for(delay - diff);
      //         // last = clock::now();
      //         update_future = std::async(std::launch::async, update_3d, last_r_cap);
      //       }
      //
      //       if (inv_idx % ((index_count - 1)/100) == 0)
      //         QMetaObject::invokeMethod(this, [=, this] {
      //           UpdateStatus(fmt::format("{:.1f} %", 100.*inv_idx/index_count));
      //         });
      //     }
      //     else {
      //       if (inv_idx % ((*pni_.connected_macro_count() - 1)/20) == 0) {
      //         // if (reference_count == pni_.connected_macro_count())
      //           if (auto diff = duration_cast<seconds>(clock::now() - last); diff < delay)
      //             std::this_thread::sleep_for(delay - diff);
      //         last = clock::now();
      //         /*update_future = */
      //         std::async(std::launch::async, update_3d, last_r_cap).wait();
      //       }
      //
      //       if (inv_idx % ((*pni_.connected_macro_count() - 1)/200) == 0)
      //         QMetaObject::invokeMethod(this, [=, this] {
      //           UpdateStatus(fmt::format("{:.1f} %", 100.*inv_idx/index_count));
      //         });
      //     }
      //
      //     
      //
      //
      //
      //
      //
      //     auto [elem, local_idx, r_cap] = queue.front();
      //
      //     queue.pop();
      //
      //     net_idx_t net_idx; // TODO
      //     if (elem == displ_elem::macro)
      //       net_idx = pni_.net(macro_idx_t(local_idx));
      //     else if (elem == displ_elem::voxel) {
      //       net_idx = pni_.net(voxel_idx_t(local_idx));
      //       darcy_invaded = true;
      //     }
      //
      //
      //     if (
      //       // et_algo::get_header(dc_props.get_entry(*net_idx)) ==
      //       // et_algo::get_header(outlet_entry)
      //       true
      //     )
      //     {
      //       // dc_context_.adjacent_edges_remove(*net_idx, dc_graph_);
      //       processed->invaded_macro_voxel[net_idx] = true;
      //
      //       last_r_cap = std::min(r_cap, last_r_cap);
      //
      //       if (elem == displ_elem::macro) {
      //         inv_volume_coefs[0] += 
      //           pn_.node_[attribs::volume][local_idx];
      //
      //         inv_volume_coefs[2] += 
      //           pn_.node_[attribs::volume][local_idx]*
      //           -props::area_of_films(theta)/props::area(pn_.node_[attribs::r_ins][local_idx]);
      //       }
      //       else if (elem == displ_elem::voxel) {
      //         inv_volume_coefs[0] += unit_darcy_pore_volume;
      //       }
      //
      //       for (auto ab : dc_graph_.edges(*net_idx))
      //         if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) { // not outlet
      //           if (!processed->explored_macro_voxel[b_net_idx]) {
      //             if (pni_.is_macro(b_net_idx)) { // macro
      //               auto b_macro_idx = pni_.macro(b_net_idx);
      //
      //               queue.insert(
      //                 displ_elem::macro, *b_macro_idx,
      //                 props::r_cap_piston_with_films(theta, pn_.node_[attribs::r_ins][*b_macro_idx]));
      //             }
      //             else { // darcy
      //               auto b_voxel_idx = pni_.voxel(b_net_idx);
      //
      //               queue.insert(displ_elem::voxel, *b_voxel_idx, darcy_r_cap_const);
      //             }
      //
      //             processed->explored_macro_voxel[b_net_idx] = true;
      //           }
      //         }
      //     }
      //   }
      //
      //   
      //
      //   QMetaObject::invokeMethod(this,
      //     [=, this] {
      //       update_3d(last_r_cap);
      //
      //       UpdateStatus(fmt::format("done {}s", duration_cast<seconds>(clock::now() - t2).count()/*/60.*/));
      //
      //       auto inv_volume = inv_volume_coefs[0] + inv_volume_coefs[2]*last_r_cap*last_r_cap;
      //
      //       add_to_plot({1 - inv_volume/total_pore_volume, 1/last_r_cap/*10000*progress*/});
      //     });
      // });
    }

    void ComputePressure() {
      using clock = std::chrono::high_resolution_clock;
      using seconds = std::chrono::seconds;

      std::cout << fmt::format("image path: {}\n\n", startup_.image.path);

      auto pnm_path = ProcessImage(startup_.image.path)/"";

      auto begin_init_time = clock::now();

      //------------------------------

      v3i processors{1};

      if (startup_.solver.decomposition)
        processors = *startup_.solver.decomposition;
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

      auto const_permeability = startup_.microporous_perm*0.001*presets::darcy_to_m2;

      
      auto cache_path = fmt::format("cache/{}-pressure-{:.2f}mD.bin",
        startup_.image.path.stem(), const_permeability/presets::darcy_to_m2*1e3);



      InitLutNodeThroat(lut_node_throat_);
      InitLutPoreSolidMicro(lut_pore_solid_micro_);

      dpl::vtk::PopulateLutRedWhiteBlue(lut_pressure_);
      dpl::vtk::PopulateLutRedWhiteBlue(lut_continuous_);

      lut_pressure_->SetTableRange(0, 1);
      // lut_pressure_->SetNanColor(0.584314, 0.552624, 0.267419, 1);

      pn_.read_from_text_file(pnm_path);

      #ifdef XPM_DEBUG_OUTPUT
        std::cout
          << fmt::format("network\n  nodes: {}\n  throats: {}\n", pn_.node_count(), pn_.throat_count())
          << (pn_.eval_inlet_outlet_connectivity() ? "(connected)" : "(disconected)") << '\n';
      #endif

      #ifdef _WIN32
        pn_.connectivity_flow_summary(startup_.solver.tolerance, startup_.solver.max_iterations);
      #else
        pn_.connectivity_flow_summary_MPI(startup.solver.tolerance, startup.solver.max_iterations);
      #endif

      std::cout << '\n';

      {
        img_.read_image(startup_.image.path, startup_.image.phases);
        img_.dim = startup_.loaded ? startup_.image.size : std::round(std::cbrt(img_.size));

        img_.read_icl_velems(pnm_path);

        auto mapped_range =
          std::ranges::subrange{img_.velem.get(), img_.velem.get() + img_.size}
        | std::views::transform([](voxel_property::velem_t x) { return *x; });

        InitLutVelems(lut_velem_, *std::ranges::max_element(mapped_range));

        img_.eval_microporous();
      }
      
      std::cout << "connectivity (isolated components)...";

      pni_.evaluate_isolated();

      std::cout << fmt::format(" done\n  macro: {}\n  voxel: {}\n\n",
        pni_.connected_macro_count(),
        pni_.connected_count() - pni_.connected_macro_count());

      using namespace presets;

      std::vector<HYPRE_Complex> pressure;
      HYPRE_Real residual = std::numeric_limits<HYPRE_Real>::quiet_NaN();
      HYPRE_Int iters = 0;


      if (startup_.use_cache && std::filesystem::exists(cache_path)) {
        std::cout << "using cached pressure\n\n";

        std::ifstream is(cache_path, std::ifstream::binary);
        is.seekg(0, std::ios_base::end);
        auto nrows = is.tellg()/sizeof(HYPRE_Complex);
        is.seekg(0, std::ios_base::beg);
        pressure.resize(nrows);
        is.read(reinterpret_cast<char*>(pressure.data()), nrows*sizeof(HYPRE_Complex));
      }
      else {
        std::cout << "decomposition...";

        auto decomposed = pni_.decompose_rows(processors);

        // for (auto i = 0; i < processors.prod(); ++i)
        //   std::cout << std::format("\nblock {}, rows {}--{}, size {}",
        //     i, decomposed.blocks[i].lower, decomposed.blocks[i].upper, decomposed.blocks[i].upper - decomposed.blocks[i].lower + 1);

        std::cout
          << " done\n\n"
          << "input matrix build...";

        auto [nrows, nvalues, input] = pni_.generate_pressure_input(decomposed, const_permeability);

        std::cout
          << " done\n\n"
          // << "pre hypre time: " << duration_cast<std::chrono::milliseconds>(clock::now() - begin_init_time).count() << "ms\n\n"
        ;

        std::cout << 
          fmt::format("save input matrix [{} MB]...", (
            nrows*(sizeof(HYPRE_Int) + sizeof(HYPRE_Complex)) +
            nvalues*(sizeof(HYPRE_BigInt) + sizeof(HYPRE_Complex)))/1024/1024);

        dpl::hypre::mpi::save(input, nrows, nvalues, decomposed.blocks, startup_.solver.tolerance, startup_.solver.max_iterations);

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

          pressure.resize(nrows);

          for (HYPRE_BigInt i = 0; i < nrows; ++i)
            pressure[decomposed.decomposed_to_net[i]] = decomposed_pressure[i];
        }

        std::cout << "pressure solved\n\n";

        if (startup_.save_cache) {
          std::filesystem::create_directory("cache");
          std::ofstream cache_stream(cache_path, std::ofstream::binary);
          cache_stream.write(reinterpret_cast<const char*>(pressure.data()), sizeof(HYPRE_Complex)*pressure.size());
          std::cout << "pressure cached\n\n";
        }
      }

      pni_.flow_summary(pressure.data(), const_permeability);
      std::cout << fmt::format("  residual: {:.4g}, iterations: {}\n\n", residual, iters);

      

      auto assembly = vtkSmartPointer<vtkAssembly>::New();

      {
        vtkSmartPointer<vtkActor> macro_actor;
        std::tie(macro_actor, macro_colors) = CreateNodeActor(pn_, lut_pressure_, 
          [&](macro_idx_t i) {
            return pni_.connected(i) ? 0/*pressure[pni_.net(i)]*/ : std::numeric_limits<double>::quiet_NaN();
          });
           
        assembly->AddPart(macro_actor);
      }

      {
        vtkSmartPointer<vtkActor> throat_actor;
        std::tie(throat_actor, throat_colors) = CreateThroatActor(pn_, lut_pressure_, [&](std::size_t i) {
          auto [l, r] = pn_.throat_[attribs::adj][i];

          return
           
            pni_.connected(l)
               ? 0/*(pressure[pni_.net(l)] + pressure[pni_.net(r)])/2*/
               : std::numeric_limits<double>::quiet_NaN();

          // return
          //   pni_.connected(l)
          //     ? pn_.inner_node(r)
          //       ? pni_.connected(r) ? 0/*(pressure[pni_.net(l)] + pressure[pni_.net(r)])/2*/ : std::numeric_limits<HYPRE_Complex>::quiet_NaN()
          //       : r == pn_.inlet() ? 1.0 : 0.0
          //     : std::numeric_limits<HYPRE_Complex>::quiet_NaN();
        });

        assembly->AddPart(throat_actor);
      }

      renderer_->AddActor(assembly);

      // std::ranges::fill(pressure, 0);
      
      // std::cout << "\n\nNetwork actor created";

      {
        auto scale_factor = /*1.0*/pn_.physical_size.x()/img_.dim.x(); // needed for vtk 8.2 floating point arithmetics
            
        img_glyph_mapper_.Init(scale_factor);
        {
          std::vector<bool> filter(img_.size);
          for (idx1d_t i = 0; i < img_.size; ++i)
            filter[i] = img_.phase[i] == microporous;

          cout << "3D faces...";

          auto t0 = clock::now();

          img_glyph_mapper_.Populate(img_.dim, pn_.physical_size/img_.dim,
            [&](idx1d_t idx1d) { return filter[idx1d]; });

          std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t0).count());
              
          dpl::sfor<6>([&](auto face_idx) {
            dpl::vtk::GlyphMapperFace<idx1d_t>& face = std::get<face_idx>(img_glyph_mapper_.faces_);

            idx1d_t i = 0;
            for (auto idx1d : face.GetIndices())
              face.GetColorArray()->SetTypedComponent(i++, 0, 
                pni_.connected(voxel_idx_t{idx1d})
                  ? /*pressure[pni_.net(v_idx)]*/ 0
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

            actor->GetProperty()->SetEdgeVisibility(/*false*/false);
            actor->GetProperty()->SetEdgeColor(v3d{0.25} /*0, 0, 0*/);
                
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



// std::filesystem::path image_path = R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Bmps-v0s255_252x252x252_6p0um.raw)";
      // parse::image_dict input_spec{
      //   .pore = 0,
      //   .solid = 1,       // dummy value, no '1' is in the image
      //   .microporous = 255, // we read actual solid as microporous
      // };

      // std::filesystem::path image_path = R"(C:\Users\dmytr\OneDrive - Imperial College London\hwu_backup\temp\images\Est-v0m2s3_500x500x500_4p0um.raw)";
      // constexpr parse::image_dict input_spec{
      //   .solid = 3,
      //   .pore = 0,
      //   .microporous = 2
      // };
     
      // HYPRE_Real tolerance = 1.e-9; HYPRE_Int max_iterations = 1000;






// auto pos = std::make_unique<v3d[]>(index_count);
        // {
        //   auto* pos_ptr = pos.get();
        //
        //   for (macro_idx i{0}; i < pn_.node_count(); ++i)
        //     if (pni_.connected(i))
        //       *pos_ptr++ = pn_.node_[attribs::pos][*i];
        //
        //   idx3d_t ijk;
        //   auto& [i, j, k] = ijk;
        //   voxel_idx idx1d;
        //
        //   auto cell_size = pn_.physical_size/img_.dim;
        //
        //   for (k = 0; k < img_.dim.z(); ++k)
        //     for (j = 0; j < img_.dim.y(); ++j)
        //       for (i = 0; i < img_.dim.x(); ++i, ++idx1d)
        //         if (pni_.connected(idx1d)) // darcy node
        //           *pos_ptr++ = cell_size*(ijk + 0.5);
        // }
        // std::sort(indices.get(), indices.get() + index_count, [&pos](idx1d_t l, idx1d_t r) { return pos[l].x() > pos[r].x(); });