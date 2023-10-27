/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * DPL is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * DPL is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <dpl/static_vector.hpp>

namespace dpl::vtk
{
  template <typename idx1d_t>
  class GlyphMapperFace
  {
  protected:
    vtkNew<vtkPoints> points_;
    std::vector<idx1d_t> original_indices_;

    vtkNew<vtkGlyph3DMapper> glyph_mapper_;
    vtkNew<vtkActor> actor_;
    vtkNew<vtkDoubleArray> color_arr_;

  public:
    auto* GetActor() const {
      return actor_.Get();
    }

    auto* GetGlyphMapper() const {
      return glyph_mapper_.Get();
    }

    auto& GetIndices() const {
      return original_indices_;
    }

    auto* GetColorArray() const {
      return color_arr_.Get();
    }
  };
  
  template<int face_idx, typename idx1d_t>
  class GlyphMapperFaceGeneric : public GlyphMapperFace<idx1d_t>
  {
    using face = face_cubic<face_idx>;
    using dims = cdims<face::dim>;
    using idx3d_t = vector_n<idx1d_t, 3>;

    static vtkSmartPointer<vtkPolyData> Quad(double half_length = 0.5) {
      auto quad = vtkSmartPointer<vtkPolyData>::New();

      vtkNew<vtkPoints> points;
      quad->SetPoints(points);
      
      vtkNew<vtkCellArray> cells;
      quad->SetPolys(cells);
        
      vector3d pos;
      pos[dims::e0] = 0;
      pos[dims::e1] = -half_length;
      pos[dims::e2] = -half_length;
      points->InsertNextPoint(pos);

      pos[dims::e1] = half_length;
      points->InsertNextPoint(pos);

      pos[dims::e2] = half_length;
      points->InsertNextPoint(pos);

      pos[dims::e1] = -half_length;
      points->InsertNextPoint(pos);

      if constexpr (face::is_upper) {
        vtkIdType indices[] = {0, 1, 2, 3};
        cells->InsertNextCell(4, indices);
      }
      else {
        vtkIdType indices[] = {3, 2, 1, 0};
        cells->InsertNextCell(4, indices);
      }

      return quad;
    }

    
  public:
    void Init(double half_length = 0.5) {
      // color_arr_->SetName("darcy_adj");

      vtkNew<vtkPolyData> polydata;
      polydata->GetPointData()->SetScalars(this->color_arr_);
      polydata->SetPoints(this->points_);
      

      this->glyph_mapper_->OrientOff();
      this->glyph_mapper_->SetScaleFactor(half_length/*1.0000*/); 
      this->glyph_mapper_->SetScaleModeToNoDataScaling();
      this->glyph_mapper_->SetInputData(polydata);
      this->glyph_mapper_->SetSourceData(Quad(/*half_length/2*/));
    }

    void Populate(const vector3i& cells, const vector3d& cell_size, const auto& filter) {
      idx3d_t map_idx{1, cells.x(), cells.x()*cells.y()};
      idx3d_t ijk;
      
      auto [e0, e1, e2] = dims::tie(ijk);
      auto [e0_count, e1_count, e2_count] = dims::tie(cells);
      
      auto adj_step = map_idx[dims::e0];

      vector3d pos;
      
      idx1d_t idx1d;
      
      for (e2 = 0; e2 < e2_count; ++e2)
        for (e1 = 0; e1 < e1_count; ++e1) {
          e0 = 0;

          idx1d = ijk.dot(map_idx);

          if constexpr (!face::is_upper) {
            if (filter(idx1d)) {
              pos[dims::e0] = 0; //(0)*cell_size[e1_dim];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];

              this->points_->InsertNextPoint(pos);
              this->original_indices_.push_back(idx1d);
            }
          }
          
          for (; e0 < e0_count - 1; ++e0) {
            idx1d = ijk.dot(map_idx);
            bool filtered = filter(idx1d);

            auto adj_idx1d = idx1d + adj_step;
            bool adj_filtered = filter(adj_idx1d);

            if (filtered != adj_filtered) {
              pos[dims::e0] = (e0 + 1.0)*cell_size[dims::e0];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];

              if constexpr (face::is_upper) {
                if (filtered) {
                  this->points_->InsertNextPoint(pos);
                  this->original_indices_.push_back(idx1d);
                }
              }
              else {
                if (adj_filtered) {
                  this->points_->InsertNextPoint(pos);
                  this->original_indices_.push_back(adj_idx1d);
                }
              }
            }
          }


          if constexpr (face::is_upper) {
            idx1d = ijk.dot(map_idx);
            
            if (filter(idx1d)) {
              pos[dims::e0] = (e0_count)*cell_size[dims::e0];
              pos[dims::e1] = (e1 + 0.5)*cell_size[dims::e1];
              pos[dims::e2] = (e2 + 0.5)*cell_size[dims::e2];
              
              this->points_->InsertNextPoint(pos);
              this->original_indices_.push_back(idx1d);
            }
          }
        }

      this->color_arr_->SetNumberOfTuples(this->points_->GetNumberOfPoints());
    }
  };


  template <typename idx1d_t>
  class ImageDataGlyphMapper
  {
  public:
    void Init(double half_length) {
      dpl::sfor<6>([=, this](auto i) {
        std::get<i>(faces_).Init(half_length);
      });
    }

    void Populate(const vector3i& dims, const vector3d& cell_size, const auto& filter) {
      dpl::psfor<6>([=, this](auto i) {
        std::get<i>(faces_).Populate(dims, cell_size, filter);
      });
    }

    std::tuple<
      GlyphMapperFaceGeneric<0, idx1d_t>,
      GlyphMapperFaceGeneric<1, idx1d_t>,
      GlyphMapperFaceGeneric<2, idx1d_t>,
      GlyphMapperFaceGeneric<3, idx1d_t>,
      GlyphMapperFaceGeneric<4, idx1d_t>,
      GlyphMapperFaceGeneric<5, idx1d_t>
    > faces_;
  };
}