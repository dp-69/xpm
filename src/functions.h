#pragma once

#include "pore_network_image.hpp"

#include <dpl/graph/et_etnte_defs.hpp>
#include <dpl/graph/dc_graph.hpp>
#include <dpl/graph/dc_context.hpp>



// #include <dpl/graph/avl_extended_augmented_tree_algorithms.hpp>
// #include <dpl/graph/cyclic_operations.hpp>
// #include <dpl/graph/dynamic_connectivity_graph.hpp>
    // #include <dpl/graph/general.hpp>



#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkCylinderSource.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3DMapper.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>

#include <algorithm>
#include <random>
#include <regex>


namespace xpm
{
  // vtkSmartPointer<vtkActor> CreateNodeActor(const pore_network_model& pnm, vtkLookupTable* lut);
  // vtkSmartPointer<vtkActor> CreateThroatActor(const pore_network_model& pnm, vtkLookupTable* lut);

  inline dpl::vector3d angles_for_j_norm(const dpl::vector3d& dir) {
    auto norm = dir.normalise();
    
    return {
      asin(norm.z())/3.141592654*180.0, 
      0.0f,                                                  
      (atan2(-norm.x(), norm.y())/3.141592654*180.0)      
    };
  }
  
  std::tuple<vtkSmartPointer<vtkActor>, vtkFloatArray*> CreateNodeActor(const pore_network& pnm, vtkLookupTable* lut, const auto& color_map) {
    using namespace attribs;

    
    vtkNew<vtkPolyData> polydata;
      
    vtkNew<vtkSphereSource> cylinder;
    cylinder->SetPhiResolution(10);
    cylinder->SetThetaResolution(10);
    cylinder->SetCenter(0, 0, 0);
    cylinder->SetRadius(1);
      
    vtkNew<vtkGlyph3DMapper> node_glyphs;
    node_glyphs->SetSourceConnection(cylinder->GetOutputPort());
    node_glyphs->SetInputData(polydata);
    node_glyphs->OrientOff();

    vtkNew<vtkFloatArray> scale_array;
    scale_array->SetName("scale");
    scale_array->SetNumberOfComponents(1);

    polydata->GetPointData()->AddArray(scale_array);
    node_glyphs->SetScaleArray(scale_array->GetName());
    node_glyphs->SetScaleModeToScaleByMagnitude();



    vtkNew<vtkFloatArray> color_array;
    color_array->SetName("color");
    color_array->SetNumberOfComponents(1);
    polydata->GetPointData()->SetScalars(color_array);

    node_glyphs->SetLookupTable(lut);
    node_glyphs->SetColorModeToMapScalars();
    node_glyphs->UseLookupTableScalarRangeOn();
    node_glyphs->SetScalarModeToUsePointData();
      
    vtkNew<vtkPoints> points;
      
    for (idx1d_t i = 0, count = pnm.node_count(); i < count; ++i) {
      points->InsertNextPoint(pnm.node_[pos][i]);
      scale_array->InsertNextTuple1(pnm.node_[r_ins][i]);
      color_array->InsertNextTuple1(color_map(macro_idx{i}));
    }
      
    polydata->SetPoints(points);
      
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(node_glyphs);

    actor->GetProperty()->SetSpecularColor(1, 1, 1);
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(75);
    actor->GetProperty()->SetAmbient(0.15);
    actor->GetProperty()->SetDiffuse(0.9);
      
    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
    return {actor, color_array};
  }

  std::tuple<vtkSmartPointer<vtkActor>, vtkFloatArray*> CreateThroatActor(const pore_network& pn, vtkLookupTable* lut, const auto& color_map) {
    using namespace attribs;

    vtkNew<vtkPolyData> polydata;
      
    vtkNew<vtkCylinderSource> cylinder;
    cylinder->SetResolution(20);
    cylinder->SetCenter(0.0, 0.5, 0.0);
    cylinder->SetRadius(1);
    cylinder->SetHeight(1);
      
    vtkNew<vtkGlyph3DMapper> throat_glyphs;
    throat_glyphs->SetSourceConnection(cylinder->GetOutputPort());
    throat_glyphs->SetInputData(polydata);
    
                     

    vtkNew<vtkFloatArray> orient_array;
    orient_array->SetName("orient");
    orient_array->SetNumberOfComponents(3);

    polydata->GetPointData()->AddArray(orient_array);
    throat_glyphs->SetOrientationArray(orient_array->GetName());
    throat_glyphs->SetOrientationModeToRotation();

    vtkNew<vtkFloatArray> scale_array;
    scale_array->SetName("scale");
    scale_array->SetNumberOfComponents(3);

    polydata->GetPointData()->AddArray(scale_array);
    throat_glyphs->SetScaleArray(scale_array->GetName());
    throat_glyphs->SetScaleModeToScaleByVectorComponents();



      
      

    vtkNew<vtkFloatArray> color_array;
    color_array->SetName("color");
    color_array->SetNumberOfComponents(1);
    polydata->GetPointData()->SetScalars(color_array);

    throat_glyphs->SetLookupTable(lut);
    throat_glyphs->SetColorModeToMapScalars();
    throat_glyphs->UseLookupTableScalarRangeOn();
    throat_glyphs->SetScalarModeToUsePointData();

    vtkNew<vtkPoints> points;

    for (std::size_t i = 0, count = pn.throat_count(); i < count; ++i)
      if (auto [l, r] = pn.throat_[adj][i];
        /*pnm.inner_node(l) && */pn.inner_node(r)) {
        auto& n0_pos = pn.node_[pos][*l];

        points->InsertNextPoint(n0_pos);
        orient_array->InsertNextTuple(angles_for_j_norm(pn.node_[pos][*r] - n0_pos));

        scale_array->InsertNextTuple(v3d{
          pn.throat_[r_ins][i],
          (pn.node_[pos][*r] - n0_pos).length(),
          pn.throat_[r_ins][i]
        });

        color_array->InsertNextTuple1(color_map(i));
      }

    polydata->SetPoints(points);
      
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(throat_glyphs);
      

    actor->GetProperty()->SetSpecularColor(1, 1, 1);
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(75);
    actor->GetProperty()->SetAmbient(0.15);
    actor->GetProperty()->SetDiffuse(0.9);

      
      
    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
    return {actor, color_array};        
  }

  


  namespace test
  {
    inline auto GenerateGraph() {
      using namespace dpl::graph;

      int vertex_count = 6;

      std::pair<int, int> edges[] = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {5, 2}, {5, 1}, {5, 0}
      };

      graph_generator gen(vertex_count);

      for (auto [l, r] : edges)
        gen.reserve(l, r);

      gen.allocate();

      for (auto [l, r] : edges)
        gen.set(l, r);

      return gen.acquire();
    }

    inline void DFS_CHECK() {
      using namespace dpl::graph;

      // for (auto i = 0; i < vertex_count; ++i) {
      //   auto [begin, end] = out_edges(g.get_vertex(i), g);
      //   for (auto iter = begin; iter != end; ++iter) {
      //     auto directed_edge = *iter;
      //   }
      // }

      auto g = GenerateGraph();

      dc_properties props{g};

      dc_context<dc_properties> context;

      context.init_with_dfs(g, props);


      // auto etnte_hdr = dc_properties::get_etnte_header(et_algo::get_header(dc_properties::get_entry(g.get_vertex(0), g)));

      #ifndef ETNTE_AS_AVL_ONLY
      for (auto etnte_entry : range<etnte_traits>(etnte_hdr))
        etnte_entry->size = 0;

      etnte_algo::populate_sizes(etnte_hdr);
      #endif

      // if (!etnte_algo::verify(etnte_hdr)) {
      //   std::cout << "Invalid recalc\n";
      // }


      auto print_clusters = [&]() {
        for (int i = 0; i < g.vertex_count(); ++i) {
          auto v = g.get_vertex(i);

          auto hdr = et_algo::get_header(dc_properties::get_entry(v, g));

          auto v_ref = dc_properties::get_vertex(
            *std::ranges::find_if(
              range<et_traits>(hdr),
              dc_properties::is_loop_edge));

          std::cout << fmt::format("vertex {} : root {}\n", props.get_idx(v), props.get_idx(v_ref));
        }
      };



      // {
      //   auto hdr = et_algo::get_header(dc_properties::get_entry(g.get_vertex(0)));
      //   print(hdr, dc_properties{g});
      //   print(dc_properties::get_etnte_header(hdr), dc_properties{g});
      //   std::cout << '\n';
      // }
      //
      // if (auto* de_for_split = g.get_directed_edge(0);
      //   dc_properties::is_tree_edge(de_for_split))
      //   etdc_context.tree_edge_split_and_reconnect(de_for_split);
      // else
      //   etdc_context.non_tree_edge_remove(de_for_split);
      //
      // {
      //   auto hdr = et_algo::get_header(dc_properties::get_entry(g.get_vertex(0)));
      //   std::cout << '\n';
      //   print(hdr, dc_properties{g});
      //   print(dc_properties::get_etnte_header(hdr), dc_properties{g});
      //   std::cout << '\n';
      // }

      print_clusters();

      for (std::size_t i = 0; i < g.vertex_count(); ++i)
        context.adjacent_edges_remove(i, g);

      std::cout << '\n';

      print_clusters();
    }

    inline void split_join_validity_check_ET_ONLY() {
      using namespace std;
      using namespace boost;
      using namespace chrono;


      using traits = dpl::graph::et_traits;
      using node = traits::node;
      using algo = dpl::graph::et_algo;
      using cyclic_op = dpl::graph::cyclic<algo>;

      using vertex = dpl::graph::vertex;
      
      random_device rand_seed;
      auto random_engine = default_random_engine(rand_seed());
      uniform_real_distribution<> unit_dist(0, 1);

      system_clock::time_point t0, t1;

      std::string text = "l lk km m mn n nj j jn nm mk ki ih hg g ge ec c cO O Oc ce ef f fe ed da a ad db b bt t tb bd d de e eg gh h hi i ik k kl";

      std::regex word_regex("(\\w+)");
      auto words_begin = std::sregex_iterator(text.begin(), text.end(), word_regex);
      auto words_end = std::sregex_iterator();

      int total_size = std::distance(words_begin, words_end); // NOLINT(cppcoreguidelines-narrowing-conversions)

      vector<node> input(total_size);

      
      

      constexpr auto display_output = false;
      constexpr auto check_validity_naive = true;

      node header_a_storage;
      auto* header_a = &header_a_storage;
      algo::init_header(header_a);


      t0 = system_clock::now();
      vector<vertex> vertices(total_size);

      for (int i = 0; i < total_size; i++) {
        dpl::graph::dc_properties::set_vertex(&input[i], &vertices[i]);
        algo::push_back(header_a, &input[i]);
      }

      t1 = system_clock::now();
      cout << boost::format("AVL creation: %i ms\n\n") % duration_cast<milliseconds>(t1 - t0).count();
      cout << "*** verification: " << algo::verify(header_a) << "***\n";


      cyclic_op::cut(header_a, &input[0]); // tip: does effectively nothing


      node headerNodeB;
      auto header_b = &headerNodeB;
      algo::init_header(header_b);


      std::size_t iter = 0;

      auto aIdx = total_size / 2;
      auto bIdx = total_size / 4;
      int cIdx;


      cout << "Binary comparison\n";

      t0 = system_clock::now();

      iter = 0;

      while (iter++ < 50000) {
        bIdx = static_cast<int>(round(unit_dist(random_engine) * (total_size - 1)));


        if (aIdx != bIdx && aIdx < bIdx != algo::less_than(&input[aIdx], &input[bIdx]))
          cout << "\nINVALID COMPARISON\n";

        aIdx = bIdx;

        if (iter % 10000 == 0) {
          t1 = system_clock::now();
          cout << boost::format("%i ms: ") % duration_cast<milliseconds>(t1 - t0).count() << iter << endl;
          t0 = t1;
        }
      }


      cout << "Split & Join\n";

      t0 = system_clock::now();

      iter = 0;


      while (iter++ < 50000)
      {
        if (total_size == 0)
          break;


        auto splitIdx = static_cast<int>(round(unit_dist(random_engine) * (total_size - 1)));


        if (iter % 10000 == 0) {
          t1 = system_clock::now();
          cout << boost::format("%i ms: ") % duration_cast<milliseconds>(t1 - t0).count() << iter << endl;
          t0 = t1;
        }

        auto splitValue = &input[splitIdx];

        auto sizeA = splitIdx;
        auto sizeB = total_size - sizeA - 1;



        algo::split_tree(header_a, splitValue, header_b);



        if (check_validity_naive) {
          if (!algo::verify(header_a))
            cout << "{A} INVALID\n";

          if (!algo::verify(header_b))
            cout << "{B} INVALID\n";
        }

        algo::join_trees(header_a, splitValue, header_b);



        if (check_validity_naive)
          if (!algo::verify(header_a))
            cout << "{A<>B} INVALID\n";
      }


      cout << "*** final verification: " << algo::verify(header_a) << "***\n";
    }

    inline void split_join_validity_check_ETNTE_ONLY() {
      using namespace std;
      using namespace boost;
      using namespace chrono;


      using traits = dpl::graph::etnte_traits;
      using node = traits::node;
      using algo = dpl::graph::etnte_algo;
      using cyclic_op = dpl::graph::cyclic<algo>;

      // using vertex = dpl::graph::vertex;
      using directed_edge = dpl::graph::directed_edge;
      
       //dpl::graph::et_cyclic_op;

      random_device rand_seed;
      // auto random_engine = default_random_engine(1500);  
      auto random_engine = default_random_engine(rand_seed());
      uniform_real_distribution<> unit_dist(0, 1);

      system_clock::time_point t0, t1;

      // auto totalSize = 2 * integral_power(10, 2 + 6);
      // auto totalSize = integral_power(10, 2 + 6);
      // auto totalSize = integral_power(10, 1 + 6);
      // auto totalSize = 9;

      //std::string text = "lh lm li ml nk jd ji gf gc cb cg fg fi fd fh ab at be ba bc ta dj df eb hl hf ij if il kn";
      //std::string text = "hg g ge ec c cO O Oc ce ef f fe ed da a ad db b bt t tb bd d de e eg gh h hi i ik k kl l lk km m mn n nj j jn nm mk ki ih";
      std::string text = "l lk km m mn n nj j jn nm mk ki ih hg g ge ec c cO O Oc ce ef f fe ed da a ad db b bt t tb bd d de e eg gh h hi i ik k kl";

      std::regex word_regex("(\\w+)");
      auto words_begin = std::sregex_iterator(text.begin(), text.end(), word_regex);
      auto words_end = std::sregex_iterator();

      int total_size = std::distance(words_begin, words_end); // NOLINT(cppcoreguidelines-narrowing-conversions)

      vector<node> input(total_size);


      constexpr auto display_output = false;
      constexpr auto check_validity_naive = true;

      // euler_tour_tree treeA;
      // auto headerA = treeA.header_ptr();
      node header_a_storage;
      auto* header_a = &header_a_storage;
      algo::init_header(header_a);


      t0 = system_clock::now();
      vector<directed_edge> edges(total_size);

      for (int i = 0; i < total_size; i++) {
        // dpl::graph::dc_properties::set_directed_edge(&input[i], &edges[i]);
        algo::push_back(header_a, &input[i]);
      }

      t1 = system_clock::now();
      cout << boost::format("AVL creation: %i ms\n\n") % duration_cast<milliseconds>(t1 - t0).count();
      // t0 = system_clock::now();
      // euler_tour_algorithms::refresh_size(headerA);
      // t1 = system_clock::now();
      // cout << boost::format("iterative: %i ms\n") % duration_cast<milliseconds>(t1 - t0).count();  
      cout << "*** verification: " << algo::verify(header_a) << "***\n";


      //     { 
      //
      //      print_latex_tree(cout, et_nt::get_parent(headerA), 0, input.data(), results.data());
      // //      for (auto& x : input) {
      // //        print_subsize(cout, &x, input.data());
      // //        cout << endl;
      // //      }
      //
      //       // et_node headerNodeB;
      //       // auto headerB = &headerNodeB;
      //       // et_algo::init_header(headerB);
      //       //
      //       //
      //       //
      //       // et_algo::split_tree(headerA, &input[26], headerB);
      //       // et_algo::erase(headerA, &input[0]);
      //       
      //
      // //      cout << "========\n\n";
      //
      //
      //       et_node headerNodeB, headerNodeC;
      //        auto headerB = &headerNodeB;
      //        auto headerC = &headerNodeC;
      //        et_algo::init_header(headerB);
      //        et_algo::init_header(headerC);
      // //
      //       et_algo::split_tree(headerA, &input[13], headerB);
      //       et_algo::split_tree(headerB, &input[39], headerC);
      //       et_algo::erase(headerA, &input[0]);
      //       et_algo::join_trees(headerC, &input[0], headerA);
      //
      // //
      // //
      // //       et_algo::split_tree(headerA, &input[7], headerB);
      // //       et_algo::split_tree(headerB, &input[24], headerC);
      // //       et_algo::erase(headerA, &input[0]);
      // //       et_algo::join_trees(headerC, &input[0], headerA);
      // //       et_algo::push_front(headerB, &input[7]);
      // //       et_algo::push_front(headerC, &input[24]);
      //
      //
      //
      //
      // //      13
      // //
      // //    39
      //
      // //      et_algo::erase(headerA, &input[0]);
      // //      et_algo::split_tree(headerA, &input[19], headerB);
      // //      cout << "========\n\n";
      // //      print_latex_tree(cout, et_nt::get_parent(headerA), 0, input.data(), results.data());
      //       cout << "========\n\n";
      //       print_latex_tree(cout, et_nt::get_parent(headerB), 0, input.data(), results.data());
      //       cout << "========\n\n";
      //       print_latex_tree(cout, et_nt::get_parent(headerC), 0, input.data(), results.data());
      //
      //      
      //
      //
      //
      //
      //
      //
      //
      //
      //
      //
      // //      auto root = et_nt::get_parent(headerA);
      // //      cout << "node {" << root->_TEXT_TO_REMOVE_ << "}\n";
      // //      print_latex_tree(cout, et_nt::get_left(root), INDENT_STEP);
      // //      print_latex_tree(cout, et_nt::get_right(root), INDENT_STEP);
      // //      cout << ';';           
      //     }

      // getchar();

      // return;


      //  eti_node int0(input[20], input[50]);
      //  eti_node int1(input[57], input[59]);
      //  eti_node int2(input[67], input[120]);
      //  eti_node int3(input[120], input[130]);
      //  eti_node int4(input[135], input[137]);
      //  eti_node int5(input[187], input[192]);
      //  
      //  eti_node inter_check(input[10], input[198]);
      //
      //  eti_tree etiTree;
      //
      //  etiTree.push_back(int0);
      //  etiTree.push_back(int1);
      //  etiTree.push_back(int2);
      //  etiTree.push_back(int3);
      //  etiTree.push_back(int4);
      //  etiTree.push_back(int5);
      //
      //
      //  auto etiHeader = etiTree.header_ptr();
      //  eti_query::align_et_to_eti(etiHeader, headerA);
      //
      //  auto found0 = eti_query::find_interval_aligned(etiHeader, input[49]);  
      //  auto vr0 = et_nt::get_vertex(eti_node_traits::get_etde0(found0))->idx_;
      //
      //  auto found1 = eti_query::find_interval_aligned(etiHeader, input[21]);
      //  auto vr1 = et_nt::get_vertex(eti_node_traits::get_etde0(found1))->idx_;
      //
      //  auto found2 = eti_query::find_interval_aligned(etiHeader, input[59]);
      //  auto found3 = eti_query::find_interval_aligned(etiHeader, input[52]);
      //  auto found4 = eti_query::find_interval_aligned(etiHeader, input[68]);
      //  auto vr4 = et_nt::get_vertex(eti_node_traits::get_etde0(found4))->idx_;
      //
      //
      //  auto found5 = eti_query::find_interval_aligned(etiHeader, input[58]);
      //  auto vr5 = et_nt::get_vertex(eti_node_traits::get_etde0(found5))->idx_;
      //
      //
      //  auto found6 = eti_query::find_interval_aligned(etiHeader, input[120]);
      //
      //  auto found7 = eti_query::find_interval_aligned(etiHeader, input[190]);
      //  auto vr7 = et_nt::get_vertex(eti_node_traits::get_etde0(found7))->idx_;
      //
      //
      //
      //  {
      //    
      //
      //    et_cyclic_op::principal_cut(headerA, eti_node_traits::get_etde0(&inter_check));
      //
      //    eti_node_ptr foundPtr;
      //
      //    while ((foundPtr = eti_query::find_overlap_aligned(etiHeader, &inter_check))) {
      //      cout << et_nt::get_vertex(eti_node_traits::get_etde0(foundPtr))->idx_ << " ";      
      //      eti_algo::erase(etiHeader, foundPtr);
      //    }
      //  }


      //  auto found2 = eti_query::find_before_insert_interval(etiTree.header_ptr(), treeA.header_ptr(), &input[19]);

      //  auto found3 = eti_query::find_before_insert_interval(etiTree.header_ptr(), treeA.header_ptr(), &input[60]);
      //  auto found5 = eti_query::find_before_insert_interval(etiTree.header_ptr(), treeA.header_ptr(), &input[55]);
      //  auto found6 = eti_query::find_before_insert_interval(etiTree.header_ptr(), treeA.header_ptr(), &input[totalSize - 1]);


      cyclic_op::cut(header_a, &input[0]); // tip: does effectively nothing


      //  t0 = system_clock::now();
      //  euler_tour_algorithms::refresh_size_recursive(treeA.header_ptr());
      //  t1 = system_clock::now();
      //  cout << format("recursive: %i ms\n") % duration_cast<milliseconds>(t1 - t0).count();  
      //  cout << "***" << euler_tour_algorithms::verify_max(treeA.header_ptr()) << "***\n";


      //  display_et(treeA);
      //  euler_tour_operations::cyclic_cut(headerA, &input[0]);
      //  display_et(treeA);


      //  vector<euler_tour_node_ptr> lowerEndStack(300);
      //  vector<euler_tour_tree> outerStack(300);
      //
      //  euler_tour_operations::subtree_remainder remainder;
      //  remainder.innerLowerEndStack = lowerEndStack.data();
      //  remainder.outerTreeStack = outerStack.data();


      //  euler_tour_tree treeB;    
      //  auto headerB = treeB.header_ptr();
      node headerNodeB;
      auto header_b = &headerNodeB;
      algo::init_header(header_b);


      size_t iter = 0;

      //  auto pseudoRandomIdx = 0;
      //
      //  auto pseudoRandomSize = 10000000;
      //  vector<int> pseudaRandomVector(pseudoRandomSize);
      //  for (auto i = 0; i < pseudoRandomSize; i++)
      //    pseudaRandomVector[i] = static_cast<int>(round(unit_dist(random_engine) * (totalSize - 1)));


      auto aIdx = total_size / 2;
      auto bIdx = total_size / 4;
      int cIdx;


      //  cout << "Ternary comparison\n";
      //
      //  t0 = system_clock::now();
      //
      //  iter = 0;
      //
      //  while (iter++ < 5000000) {
      ////    if (pseudoRandomIdx == pseudoRandomSize)
      ////      pseudoRandomIdx = 0;
      //
      //    cIdx = static_cast<int>(round(unit_dist(random_engine) * (totalSize - 1)));
      //
      //    euler_tour_interval_node_traits::less_than_low_high(&input[aIdx], &input[bIdx], &input[cIdx]);
      //
      //
      //    aIdx = bIdx;
      //    bIdx = cIdx;
      //        
      //
      //    if (iter%1000000 == 0) {      
      //      t1 = system_clock::now();      
      //      cout << format("%i ms: ") % duration_cast<milliseconds>(t1 - t0).count() << iter << endl;
      //      t0 = t1;     
      //    }
      //
      //   
      //  }


      cout << "Binary comparison\n";

      t0 = system_clock::now();

      iter = 0;

      while (iter++ < 50000) {
        bIdx = static_cast<int>(round(unit_dist(random_engine) * (total_size - 1)));


        if (aIdx != bIdx && aIdx < bIdx != algo::less_than(&input[aIdx], &input[bIdx]))
          cout << "\nINVALID COMPARISON\n";

        aIdx = bIdx;

        if (iter % 10000 == 0) {
          t1 = system_clock::now();
          cout << boost::format("%i ms: ") % duration_cast<milliseconds>(t1 - t0).count() << iter << endl;
          t0 = t1;
        }
      }


      cout << "Split & Join\n";

      t0 = system_clock::now();

      iter = 0;


      while (iter++ < 50000)
      //  while (iter++ >= 0)       
      {
        //    if (pseudoRandomIdx == pseudoRandomSize)
        //      pseudoRandomIdx = 0;


        if (total_size == 0)
          break;


        auto splitIdx = static_cast<int>(round(unit_dist(random_engine) * (total_size - 1)));


        //    auto splitIdx = 0;

        //    auto splitIdx = pseudaRandomVector[pseudoRandomIdx++];                         


        //    auto splitIdx = 1195710;


        //    if (splitIdx == 0 || splitIdx == totalSize - 1)
        //      splitIdx = 5;


        //    cout << endl << format("total: %i; splitIdx: %i") % totalSize % splitIdx << endl;
        if (iter % 10000 == 0) {
          //      if (!avl_augmented_tree_algorithms::verify_max(treeA.header_ptr()))
          //        cout << "{A} INVALID\n";

          //      if (!avl_augmented_tree_algorithms::verify_max(treeB.header_ptr()))
          //        cout << "{B} INVALID\n";

          t1 = system_clock::now();
          cout << boost::format("%i ms: ") % duration_cast<milliseconds>(t1 - t0).count() << iter << endl;
          t0 = t1;
        }

        //    if (iter == 6063) {
        //      int www  = 3;
        //    }

        // if (display_output) {
        //   cout << "\nTREE {A}\n";
        //   //      for (auto& x : treeA)
        //   ////        cout << value_traits::to_node_ptr(x)->MY_V0 << " ";        
        //   //        cout << x.v0 << " ";
        //   cout << endl;
        // }


        auto splitValue = &input[splitIdx];

        auto sizeA = splitIdx;
        auto sizeB = total_size - sizeA - 1;


        // if (check_validity_naive)
        //   if (!et_algo::verify(header_a))
        //     cout << "{A} INVALID\n";
        //
        // if (check_validity_naive)
        //   if (!et_algo::verify(headerB))
        //     cout << "{B} INVALID\n";


        algo::split_tree(header_a, splitValue, header_b);


        //  //    auto rootASize = node_traits::get_size(treeA.root().pointed_node());
        //      if (sizeA != euler_tour_node_traits::get_size(et_nt::get_parent(headerA)))
        //      {
        //  //      cout << splitIdx;
        //  //      getchar();
        //        cout << "{A} INVALID SIZE\n";        
        //      }
        //
        //  //    auto rootBSize = node_traits::get_size(treeB.root().pointed_node());
        //      if (sizeB != euler_tour_node_traits::get_size(et_nt::get_parent(headerB)))
        //        cout << "{B} INVALID SIZE\n";


        if (check_validity_naive) {
          if (auto* root = traits::get_parent(header_a); !(
            algo::verify(header_a)

            #ifndef ETNTE_AS_AVL_ONLY
              && algo::calculate_subtree_size(root) == (root ? traits::get_size(root) : 0)
            #endif
            ))

            cout << "{A} INVALID\n";

          if (auto* root = traits::get_parent(header_b); !(
            algo::verify(header_b)

            #ifndef ETNTE_AS_AVL_ONLY
              && algo::calculate_subtree_size(root) == (root ? traits::get_size(root) : 0)
            #endif
            ))
            cout << "{B} INVALID\n";

          //      auto calcSizeA = treeA.size();            

          //      if (calcSizeA != sizeA || calcSizeA != euler_tour_node_traits::get_size(treeA.root().pointed_node()))
          //        cout << "{A} INVALID SIZE\n";    

          //      auto calcSizeB = treeB.size();

          //      if (calcSizeB != sizeB || calcSizeB != euler_tour_node_traits::get_size(treeB.root().pointed_node()))
          //        cout << "{B} INVALID SIZE\n";
        }

        // if (display_output) {
        //   cout << endl << "\nSPLIT\n";
        //
        //   cout << "TREE {A}\n";
        //   //      for (auto& x : treeA)
        //   ////        cout << value_traits::to_node_ptr(x)->MY_V0 << " ";
        //   //        cout << x.v0 << " ";
        //
        //   //      cout << "\n\n___" << splitValue.v0 << "___\n\n";
        //   //      cout << "\n\n___" << value_traits::to_node_ptr(splitValue)->MY_V0 << "___\n\n";
        //
        //   cout << "TREE {B}\n";
        //   //      for (auto& x : treeB)
        //   ////        cout << value_traits::to_node_ptr(x)->MY_V0 << " ";
        //   //      cout << x.v0 << " ";
        // }


        // if (check_validity_naive) {
        //   auto preHeightA = et_algo::node_height(et_traits::get_parent(header_a));
        //   auto preHeightB = et_algo::node_height(et_traits::get_parent(headerB));
        //
        //   // auto heightIncrQ = et_algo::join_trees(header_a, splitValue, headerB);
        //
        //   auto diff =
        //     et_algo::node_height(et_traits::get_parent(header_a)) -
        //     max(preHeightA, preHeightB);
        //
        //   // if (diff > 1 || heightIncrQ && diff == 0 || !heightIncrQ && diff == 1)
        //   //   cout << "MERGE HEIGHT ERROR\n";
        //
        //   if (!et_algo::verify(header_a))
        //     cout << "{A<>B} INVALID\n";
        //
        //   //      auto treeAsize = treeA.size();
        //   //      auto rootASize = euler_tour_node_traits::get_size(treeA.root().pointed_node());
        //
        //   //      if (treeAsize != totalSize || totalSize != rootASize)
        //   //        cout << "{A<>B} INVALID SIZE\n";
        //
        //   //      auto calcSizeB = treeB.size();
        //
        //   //      if (calcSizeB != 0)
        //   //        cout << "{B} INVALID SIZE\n";
        // }

        algo::join_trees(header_a, splitValue, header_b);

        //      if (totalSize != euler_tour_node_traits::get_size(euler_tour_node_traits::get_parent(headerA)))
        //        cout << "{A<>B} INVALID SIZE\n";       


        if (check_validity_naive)
          if (!algo::verify(header_a))
            cout << "{A<>B} INVALID\n";

        // if (display_output) {
        //   cout << "\nMERGED TREE {A<>B}\n";
        //
        //   //      for (auto& x : treeA)                
        //   //        cout << x.v0 << " ";
        //   ////      cout << value_traits::to_node_ptr(splitValue)->MY_V0 << " ";
        //
        //   cout << endl;
        // }


        //    auto valid = euler_tour_directed_edge_algorithms::verify(headerA);

        //    euler_tour_directed_edge_algorithms::erase(headerA, splitValue);

        //    myPool.destroy(splitValue);

        //    auto valid2 = euler_tour_directed_edge_algorithms::verify(headerA);

        //    --totalSize;
        //    for (auto i = splitIdx; i < totalSize; i++)
        //      input[i] = input[i + 1];
      }


      cout << "*** final verification: " << algo::verify(header_a) << "***\n";
    }
  }
}



