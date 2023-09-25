#pragma once

#include "pore_network_image.hpp"

namespace xpm {
  class invasion_task {
    pore_network_image* pni_;
    pore_network* pn_;
    image_data* img_;

    dpl::graph::dc_graph dc_graph_;
    std::unordered_map<dpl::graph::dc_graph::edge_descriptor, std::size_t> de_to_throat_;

    dpl::graph::dc_context<dpl::graph::dc_properties> dc_context_;

    dpl::strong_array<net_tag, bool> invaded_macro_voxel_;
    std::vector<bool> invaded_throat_;

    std::vector<dpl::vector2d> pc_curve_;

    idx1d_t inv_idx_ = 0;

    double last_r_cap_ = std::numeric_limits<double>::max();

    bool darcy_invaded_ = false;
    bool finished_ = false;


    void add_to_pc_curve(const dpl::vector2d& p) {
      if (pc_curve_.size() > 1
        && std::abs(pc_curve_[pc_curve_.size() - 1].y() - p.y()) < 1e-6
        && std::abs(pc_curve_[pc_curve_.size() - 2].y() - p.y()) < 1e-6)
        pc_curve_[pc_curve_.size() - 1].x() = p.x();
      else
        pc_curve_.push_back(p);
    }

  public:
    explicit invasion_task(pore_network_image& pni)
      : pni_{&pni}, pn_{&pni.pn()}, img_{&pni.img()} {}

    auto& graph() {
      return dc_graph_;
    }

    bool invaded(voxel_idx_t i) {
      return invaded_macro_voxel_[pni_->net(i)];
    }

    bool invaded(macro_idx_t i) {
      return invaded_macro_voxel_[pni_->net(i)];
    }

    bool invaded(std::size_t i) {
      return invaded_throat_[i];
    }

    auto& pc_curve() {
      return pc_curve_;
    }

    auto last_r_cap() {
      return last_r_cap_;
    }

    auto inv_idx() const {
      return inv_idx_;
    }

    bool finished() const {
      return finished_;
    }

    void init() {
      using clock = std::chrono::high_resolution_clock;
      using seconds = std::chrono::seconds;

      auto t0 = clock::now();
      std::cout << "graph...";

      std::tie(dc_graph_, de_to_throat_) = pni_->generate_dc_graph<true>();
      std::cout << fmt::format(" done {}s\n  {:L} vertices\n  {:L} edges\n\n",
        duration_cast<seconds>(clock::now() - t0).count(),
        dc_graph_.vertex_count(),
        dc_graph_.edge_count());

      auto t1 = clock::now();
      std::cout << "euler tour...";
      dc_context_.init_with_dfs(dc_graph_, dpl::graph::dc_properties{dc_graph_});
      std::cout << fmt::format(" done {}s\n\n", duration_cast<seconds>(clock::now() - t1).count());

      std::cout << "full decremental connectivity... (async)\n";

      pc_curve_.reserve(1000);
    }

    void launch(double darcy_porosity, double theta, std::span<const dpl::vector2d> darcy_pc_to_sw) {
      auto index_count = *pni_->connected_count();

      using props = hydraulic_properties::equilateral_triangle_properties;

      invaded_macro_voxel_.resize(pni_->connected_count());
      invaded_throat_.resize(pn_->throat_count());

      


      double darcy_r_cap_const;

      if (darcy_pc_to_sw.empty()) {
        auto min_r_cap_throat = std::numeric_limits<double>::max();
      
        for (std::size_t i{0}; i < pn_->throat_count(); ++i)
          if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && pni_->connected(l))
            min_r_cap_throat = std::min(min_r_cap_throat, pn_->throat_[attribs::r_ins][i]);
        min_r_cap_throat = 0.95*props::r_cap_piston_with_films(theta, min_r_cap_throat);
      
        darcy_r_cap_const = min_r_cap_throat; // TODO: USE DATA FROM THE CURVE!
      }
      else
        darcy_r_cap_const = 1/darcy_pc_to_sw.front().x();



      dpl::strong_array<net_tag, bool> explored(pni_->connected_count());
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue queue;

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inlet() == r) // inlet macro nodes TODO: voxels inlet needed
          if (pni_->connected(l)) {
            queue.insert(displ_elem::throat, i, props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][i]));
            explored_throat[i] = true;
          }

      net_idx_t outlet_idx{index_count};

      dpl::graph::dc_properties dc_props{dc_graph_};
      dpl::graph::et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);

      dpl::vector3d inv_volume_coefs{0};
      idx1d_t inv_darcy_count{0};

      last_r_cap_ = queue.front().radius_cap;

      auto total_pore_volume = pni_->eval_total_pore_volume(darcy_porosity);

      dpl::vector2d last_pc_point{2, -1};

      auto unit_darcy_pore_volume = (pn_->physical_size/img_->dim).prod()*darcy_porosity;

      

      auto eval_inv_volume = [&] {
        return
          inv_volume_coefs[0] + inv_volume_coefs[2]*last_r_cap_*last_r_cap_
        + (1 - (darcy_pc_to_sw.empty() ? 0.0 : solve(darcy_pc_to_sw, 1/last_r_cap_, dpl::extrapolant::flat)))*unit_darcy_pore_volume*inv_darcy_count
        ;
      };

      for (inv_idx_ = 0; !queue.empty(); ++inv_idx_) {
        // if (inv_idx_ < 50)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{100});
        // else if (inv_idx_ < 100)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{50});
        // else if (inv_idx_ < 1000)
        //   std::this_thread::sleep_for(std::chrono::milliseconds{10});

        auto inv_volume = eval_inv_volume();

        if (dpl::vector2d pc_point{1 - inv_volume/total_pore_volume, 1/last_r_cap_};
          std::abs(last_pc_point.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point.y() - pc_point.y()) > 1e4) {
          last_pc_point = pc_point;
          add_to_pc_curve(pc_point);
        }


        auto [elem, local_idx, r_cap] = queue.front();

        queue.pop();



        if (elem == displ_elem::macro) {
          auto net_idx = pni_->net(macro_idx_t(local_idx));
          invaded_macro_voxel_[net_idx] = true;

          last_r_cap_ = std::min(r_cap, last_r_cap_);

          inv_volume_coefs[0] += pn_->node_[attribs::volume][local_idx];

          inv_volume_coefs[2] += pn_->node_[attribs::volume][local_idx]*
            -props::area_of_films(theta)/props::area(pn_->node_[attribs::r_ins][local_idx]);

          for (auto ab : dc_graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) // not outlet
              if (pni_->is_macro(b_net_idx)) { // macro (throat)
                if (auto throat_idx = de_to_throat_[ab]; !explored_throat[throat_idx]) {
                  queue.insert(displ_elem::throat, throat_idx, props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]));
                  explored_throat[throat_idx] = true;
                }
              }
              else { // darcy
                if (!explored[b_net_idx]) {
                  auto b_voxel_idx = pni_->voxel(b_net_idx);
                  queue.insert(displ_elem::voxel, *b_voxel_idx, /*1/microprs_pc.back().y()*/ darcy_r_cap_const);
                  explored[b_net_idx] = true;
                }
              }
        }
        else if (elem == displ_elem::voxel) {
          darcy_invaded_ = true;

          auto net_idx = pni_->net(voxel_idx_t(local_idx));
          invaded_macro_voxel_[net_idx] = true;

          last_r_cap_ = std::min(r_cap, last_r_cap_);

          ++inv_darcy_count;

          for (auto ab : dc_graph_.edges(*net_idx))
            if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) // not outlet
              if (!explored[b_net_idx])
                if (pni_->is_macro(b_net_idx)) { // macro
                  auto b_macro_idx = pni_->macro(b_net_idx);
                  queue.insert(displ_elem::macro, *b_macro_idx, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*b_macro_idx]));
                  explored[b_net_idx] = true;
                }
                else { // darcy
                  auto b_voxel_idx = pni_->voxel(b_net_idx);
                  queue.insert(displ_elem::voxel, *b_voxel_idx, /*1/microprs_pc.back().y()*/ darcy_r_cap_const);
                  explored[b_net_idx] = true;
                }
        }
        else { // throat
          invaded_throat_[local_idx] = true; // TODO

          auto [l, r] = pn_->throat_[attribs::adj][local_idx]; 

          if (pn_->inner_node(r)/* && pni_->connected(l)*/) {
            last_r_cap_ = std::min(r_cap, last_r_cap_);

            inv_volume_coefs[0] += 
              pn_->throat_[attribs::volume][local_idx];

            inv_volume_coefs[2] += 
              pn_->throat_[attribs::volume][local_idx]*
              -props::area_of_films(theta)/props::area(pn_->throat_[attribs::r_ins][local_idx]);


            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(displ_elem::macro, *r, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*r]));
              explored[r_net] = true;
            }
          }

          if (auto l_net = pni_->net(l); !explored[l_net]) {
            queue.insert(displ_elem::macro, *l, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*l]));
            explored[l_net] = true;
          }
        }
      }

      for (auto mult : {1.004, 1.008, 1.016, 1.032, 1.064, 1.128}) {
        // std::this_thread::sleep_for(std::chrono::milliseconds{1250});

        auto inv_volume = eval_inv_volume();
        add_to_pc_curve({1 - inv_volume/total_pore_volume, 1/last_r_cap_});
        last_r_cap_ /= mult;
      }

      // for (auto i = 0; i < 5; ++i) {
      //   auto inv_volume = eval_inv_volume();
      //   add_to_plot({1 - inv_volume/total_pore_volume, 1/last_r_cap});
      //   last_r_cap *= 0.925;
      // }

      finished_ = true;
    }
  };
}