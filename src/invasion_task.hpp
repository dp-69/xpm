#pragma once

#include "pore_network_image.hpp"

namespace xpm {
  class invasion_task {
    pore_network_image* pni_;
    pore_network* pn_;

    dpl::graph::dc_graph dc_graph_;
    std::unordered_map<dpl::graph::dc_graph::edge_descriptor, std::size_t> de_to_throat_;

    dpl::graph::dc_context<dpl::graph::dc_properties> dc_context_;

    bool darcy_invaded_ = false;

    dpl::strong_array<net_tag, bool> invaded_macro_voxel_;
    std::vector<bool> invaded_throat_;

    std::vector<dpl::vector2d> pc_curve_;

    idx1d_t inv_idx_ = 0;

    bool finished_ = false;

  public:
    explicit invasion_task(pore_network_image& pni)
      : pni_{&pni}, pn_{&pni.pn()} {}

    auto& graph() {
      return dc_graph_;
    }

    auto& invaded_array() {
      return invaded_macro_voxel_;
    }

    auto& invaded_throat() {
      return invaded_throat_;
    }

    auto& pc_curve() {
      return pc_curve_;
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

    void launch(double microprs_porosity, std::span<dpl::vector2d> microprs_pc) {
      auto theta = 0.0;

      auto index_count = *pni_->connected_count();

      using props = hydraulic_properties::equilateral_triangle_properties;

      invaded_macro_voxel_.resize(pni_->connected_count());
      invaded_throat_.resize(pn_->throat_count());

      double darcy_r_cap_const = 
        props::r_cap_piston_with_films(theta, std::ranges::min(pn_->node_.range(attribs::r_ins)))*0.95;

      std::cout << fmt::format("max macro Pc: {}\n", 1/darcy_r_cap_const);

      auto min_r_cap_throat = std::numeric_limits<double>::max();

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && pni_->connected(l))
          min_r_cap_throat = std::min(min_r_cap_throat, pn_->throat_[attribs::r_ins][i]);

      min_r_cap_throat = 0.95*props::r_cap_piston_with_films(theta, min_r_cap_throat);

      std::cout << fmt::format("max throat Pc: {}\n", 1/min_r_cap_throat);

      darcy_r_cap_const = min_r_cap_throat;


      dpl::strong_array<net_tag, bool> explored{pni_->connected_count()};
      std::vector<bool> explored_throat(pn_->throat_count());

      displ_queue queue;

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inlet() == r) // inlet macro nodes TODO: voxels inlet needed
          if (pni_->connected(l)) {

            queue.insert(displ_elem::throat, i, props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][i]));
            explored_throat[i] = true;

            // std::cout << fmt::format("inlet throat volume: {}\n", pn_->throat_[attribs::volume][i]);

            // queue.insert(displ_elem::macro, *l, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*l]));
            // explored[pni_->net(l)] = true;
          }

      dpl::graph::dc_properties dc_props{dc_graph_};

      net_idx_t outlet_idx{index_count};
      dpl::graph::et_traits::node_ptr outlet_entry = dc_props.get_entry(*outlet_idx);

      v3d inv_volume_coefs{0};
      idx1d_t inv_darcy_count{0};

      auto last_r_cap = queue.front().radius_cap;

      auto total_pore_volume = 0.0;

      for (macro_idx_t i{0}; i < pn_->node_count(); ++i)
        if (pni_->connected(i))
          total_pore_volume += pn_->node_[attribs::volume][*i];

      for (std::size_t i{0}; i < pn_->throat_count(); ++i)
        if (auto [l, r] = pn_->throat_[attribs::adj][i]; pn_->inner_node(r) && pni_->connected(l))
          total_pore_volume += pn_->throat_[attribs::volume][i];
        
      auto unit_darcy_pore_volume = (pn_->physical_size/pni_->img().dim).prod()*microprs_porosity;

      for (voxel_idx_t i{0}; i < pni_->img().size; ++i)
        if (pni_->connected(i))
          total_pore_volume += unit_darcy_pore_volume;

      dpl::vector2d last_pc_point{2, -1};


      auto add_to_plot = [this](dpl::vector2d p) {
        if (pc_curve_.size() > 1
          && std::abs(pc_curve_[pc_curve_.size() - 1].y() - p.y()) < 1e-6
          && std::abs(pc_curve_[pc_curve_.size() - 2].y() - p.y()) < 1e-6)
          pc_curve_[pc_curve_.size() - 1].x() = p.x();
        else
          pc_curve_.push_back(p);
      };


      using namespace std::ranges::views;
      auto query = microprs_pc | reverse | transform([](dpl::vector2d p) { return dpl::vector2d{p.y(), p.x()}; });
      const std::vector<dpl::vector2d> pc_inverse{query.begin(), query.end()};

      auto eval_inv_volume = [&]() {
        return
          inv_volume_coefs[0] + inv_volume_coefs[2]*last_r_cap*last_r_cap
        // + (1 - solve(std::span{pc_inverse}, 1/last_r_cap))*unit_darcy_pore_volume*inv_darcy_count
        ;
      };

      for (inv_idx_ = 0; !queue.empty(); ++inv_idx_) {
        // std::this_thread::sleep_for(std::chrono::milliseconds{1});

        auto inv_volume = eval_inv_volume();

        if (dpl::vector2d pc_point{1 - inv_volume/total_pore_volume, 1/last_r_cap};
          std::abs(last_pc_point.x() - pc_point.x()) > 0.05 || std::abs(last_pc_point.y() - pc_point.y()) > 1e4) {
          last_pc_point = pc_point;
          add_to_plot(pc_point);
        }


        auto [elem, local_idx, r_cap] = queue.front();

        queue.pop();

        net_idx_t net_idx; 
        if (elem == displ_elem::macro)
          net_idx = pni_->net(macro_idx_t(local_idx));
        else if (elem == displ_elem::voxel) {
          net_idx = pni_->net(voxel_idx_t(local_idx));
          darcy_invaded_ = true;
        }
        else if (elem == displ_elem::throat) {
        }


        if (
          // dpl::graph::et_algo::get_header(dc_props.get_entry(*net_idx)) ==
          // dpl::graph::et_algo::get_header(outlet_entry)
          true
        )
        {
          // dc_context_.adjacent_edges_remove(*net_idx, dc_graph_);

          if (elem == displ_elem::macro || elem == displ_elem::voxel)
            invaded_macro_voxel_[net_idx] = true;

          if (elem == displ_elem::macro) {
            last_r_cap = std::min(r_cap, last_r_cap);

            inv_volume_coefs[0] += 
              pn_->node_[attribs::volume][local_idx];

            inv_volume_coefs[2] += 
              pn_->node_[attribs::volume][local_idx]*
              -props::area_of_films(theta)/props::area(pn_->node_[attribs::r_ins][local_idx]);
          }
          else if (elem == displ_elem::voxel) {
            last_r_cap = std::min(r_cap, last_r_cap);

            // inv_volume_coefs[0] += unit_darcy_pore_volume;
            ++inv_darcy_count;
          }
          else if (elem == displ_elem::throat) {
            if (auto [l, r] = pn_->throat_[attribs::adj][local_idx]; pn_->inner_node(r)/* && pni_->connected(l)*/) {
              invaded_throat_[local_idx] = true;

              last_r_cap = std::min(r_cap, last_r_cap);

              inv_volume_coefs[0] += 
                pn_->throat_[attribs::volume][local_idx];

              inv_volume_coefs[2] += 
                pn_->throat_[attribs::volume][local_idx]*
                -props::area_of_films(theta)/props::area(pn_->throat_[attribs::r_ins][local_idx]);
            }
          }





          if (elem == displ_elem::throat) {
            auto [l, r] = pn_->throat_[attribs::adj][local_idx];

            if (auto l_net = pni_->net(l); !explored[l_net]) {
              queue.insert(
                displ_elem::macro, *l, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*l]));
              explored[l_net] = true;
            }

            if (auto r_net = pni_->net(r); !explored[r_net]) {
              queue.insert(
                displ_elem::macro, *r, props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*r]));
              explored[r_net] = true;
            }
          }
          else {
            for (auto ab : dc_graph_.edges(*net_idx))
              if (net_idx_t b_net_idx{target(ab, dc_graph_)}; b_net_idx != outlet_idx) { // not outlet
                if (!explored[b_net_idx]) {
                  if (pni_->is_macro(b_net_idx)) { // macro
                    if (pni_->is_macro(net_idx)) {
                      auto throat_idx = de_to_throat_[ab];



                      auto pc = 1/props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]);

                      auto [a, b] = pn_->throat_[attribs::adj][throat_idx];

                      if (props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]) < darcy_r_cap_const)
                        std::cout << fmt::format("error value: {}, [a, b]: {}, {}\n", pc, a, b);

                      queue.insert(displ_elem::throat, throat_idx,
                        props::r_cap_piston_with_films(theta, pn_->throat_[attribs::r_ins][throat_idx]));

                      explored_throat[throat_idx] = true;
                    }
                    else {
                      auto b_macro_idx = pni_->macro(b_net_idx);

                      queue.insert(
                        displ_elem::macro, *b_macro_idx,
                        props::r_cap_piston_with_films(theta, pn_->node_[attribs::r_ins][*b_macro_idx]));
                      
                      explored[b_net_idx] = true;
                    }
                  }
                  else { // darcy
                    auto b_voxel_idx = pni_->voxel(b_net_idx);
                    queue.insert(displ_elem::voxel, *b_voxel_idx, /*1/microprs_pc.back().y()*/  darcy_r_cap_const);

                    explored[b_net_idx] = true;
                  }
                }
              }
          }
        }
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