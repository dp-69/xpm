/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <dpl/static_vector.hpp>

#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCommand.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkViewport.h>

#ifdef __cpp_lib_format
  #include <format>
#else
  #include <fmt/format.h>
#endif

#include <array>
#include <numbers>

namespace dpl::vtk
{
  template <int face>
  class FaceProperties
  {
    friend class TidyAxes;
    
    static constexpr sface_t<face> sface{};
    static constexpr cdims rel{sface};
    
    vtkNew<vtkPoints> points_;
    vtkNew<vtkCellArray> lines_;
    vtkNew<vtkPolyData> poly_data_;
    vtkNew<vtkPolyDataMapper> mapper_;

    vtkRenderer* viewport_;
    
    vector3d color_;

    const double* bounds_;
  
    vtkNew<vtkActor> gridlines_actor_;

    void Init(vtkRenderer* renderer, vector3d color, float line_width) {
      poly_data_->SetPoints(points_);
      poly_data_->SetLines(lines_);
     
      mapper_->SetInputData(poly_data_);

      gridlines_actor_->SetMapper(mapper_);
      gridlines_actor_->UseBoundsOff();
      gridlines_actor_->PickableOff();
      
      gridlines_actor_->GetProperty()->SetColor(color);
      gridlines_actor_->GetProperty()->SetLineWidth(line_width);

      renderer->AddActor(gridlines_actor_);
      viewport_ = renderer;
    }

    template <int span_dim_idx>
    void BuildGridlines(sdim<3, span_dim_idx>, double min, double step, int count) {
      static constexpr auto run_dim_idx = third<rel.i, span_dim_idx>();  // NOLINT(CppInconsistentNaming)
      
      vector3d r;

      r[rel.i] = bounds_[sface];

      auto point_idx = points_->GetNumberOfPoints();

      for (int i = 0; i < count; ++i) {
        r[run_dim_idx] = min + step*i;
        
        r[span_dim_idx] = bounds_[span_dim_idx*2];  // NOLINT(bugprone-implicit-widening-of-multiplication-result)
        points_->InsertNextPoint(r);

        r[span_dim_idx] = bounds_[span_dim_idx*2 + 1];
        points_->InsertNextPoint(r);

        lines_->InsertNextCell(2);
        lines_->InsertCellPoint(point_idx++);
        lines_->InsertCellPoint(point_idx++);
      }
    }

    
    void BuildGridlines(vector3d min, vector3d step, vector3i count) {
      points_->Initialize();
      lines_->Initialize();

      BuildGridlines(rel.k, min[rel.j], step[rel.j], count[rel.j]);
      BuildGridlines(rel.j, min[rel.k], step[rel.k], count[rel.k]);

      mapper_->Update();
    }

    void RefreshVisibility() {
      const vector2i_map viewport_size{viewport_->GetSize()};
      
      if (viewport_size.sum() == 0)
        return;

      if (auto* camera = viewport_->GetActiveCamera();
        (vector3d_map{camera->GetFocalPoint()} - vector3d_map{camera->GetPosition()})
        .normalise().dot(sface.normal) > 0.2) { // 0.2 is a manual coefficient, a more robust solution is preferred
        gridlines_actor_->SetVisibility(false);
        return;
      }
      
      vector2d face_triangle_display[3];
      
      vtkNew<vtkCoordinate> coord;
      coord->SetCoordinateSystemToWorld();
      
      vector3d r;

      r[rel.i] = bounds_[sface];
      r[rel.j] = bounds_[rel.j*2];
      r[rel.k] = bounds_[rel.k*2];
      coord->SetValue(r);
      face_triangle_display[0] = coord->GetComputedDoubleViewportValue(viewport_);

      r[rel.k] = bounds_[rel.k*2 + 1];
      coord->SetValue(r);
      face_triangle_display[1] = coord->GetComputedDoubleViewportValue(viewport_);
      
      r[rel.j] = bounds_[rel.j*2 + 1];
      coord->SetValue(r);
      face_triangle_display[2] = coord->GetComputedDoubleViewportValue(viewport_);

      
      gridlines_actor_->SetVisibility(
        sface.non_zero_component*
        operation::cross::z(
          face_triangle_display[1] - face_triangle_display[0],
          face_triangle_display[2] - face_triangle_display[1]) < 0);
    }
  };

    
  class TidyAxes
  {
    bool visibility_ = true;

    std::array<double, 6> bounds_;

    std::tuple<
      FaceProperties<0>,
      FaceProperties<1>,
      FaceProperties<2>,
      FaceProperties<3>,
      FaceProperties<4>,
      FaceProperties<5>> face_properties_;

    vtkRenderer* renderer_;
    vector3d scale_{1};

    struct {
      vtkNew<vtkPoints> points_;
      vtkNew<vtkCellArray> lines_;
      vtkNew<vtkPolyData> poly_data_;
      vtkNew<vtkPolyDataMapper> mapper_;
      vtkNew<vtkUnsignedCharArray> colors_;
      vtkNew<vtkActor> actor_;
    } axes_box_;

    struct {
      vtkNew<vtkPoints> points_;
      vtkNew<vtkCellArray> lines_;
      vtkNew<vtkPolyData> poly_data_;
      vtkNew<vtkPolyDataMapper2D> mapper_;
      // vtkNew<vtkUnsignedCharArray> colors_;
      vtkNew<vtkActor2D> actor_;
    } axes_ticks_;

    vector3d min_;
    vector3d step_;
    vector3i count_;
    
    static constexpr auto PI = std::numbers::pi_v<double>;

    static vector2i GetAlignment(vector2d v) {
      static constexpr auto vertical_narrow_factor = 1.2*PI/16;//PI/16;
      
      auto angle = std::atan2(-v.y(), -v.x()) + PI;
      
      if (angle < PI/8) // Left
        return {VTK_TEXT_LEFT, VTK_TEXT_CENTERED};
      
      if (angle < 3*PI/8 + vertical_narrow_factor) // Bottom Left
        return {VTK_TEXT_LEFT, VTK_TEXT_BOTTOM};
      
      if (angle < 5*PI/8 - vertical_narrow_factor) // Bottom
        return {VTK_TEXT_CENTERED, VTK_TEXT_BOTTOM};
      
      if (angle < 7*PI/8) // Bottom Right
        return {VTK_TEXT_RIGHT, VTK_TEXT_BOTTOM};
      
      if (angle < 9*PI/8) // Right
        return {VTK_TEXT_RIGHT, VTK_TEXT_CENTERED};
      
      if (angle < 11*PI/8 + vertical_narrow_factor) // Top Right
        return {VTK_TEXT_RIGHT, VTK_TEXT_TOP};
      
      if (angle < 13*PI/8 - vertical_narrow_factor) // Top
        return {VTK_TEXT_CENTERED, VTK_TEXT_TOP};
      
      if (angle < 15*PI/8) // Top Left
        return {VTK_TEXT_LEFT, VTK_TEXT_TOP};

      // Right
      return {VTK_TEXT_LEFT, VTK_TEXT_CENTERED};
    }

    template <int face0, int face1>
    void RefreshAxis(sface_t<face0> f0 = {}, sface_t<face1> f1 = {}) {
      if (!Face(f0).gridlines_actor_->GetVisibility() &&
          !Face(f1).gridlines_actor_->GetVisibility())
        return;

      // #ifdef __cpp_lib_format
      //   auto format = "{:" + format_ + "}";
      // #else
      //   // std::string format = "%" + format_;
      //   // const char* format_ptr = format.c_str();
      // #endif
      auto format = "{:" + format_ + "}";

      
      static constexpr auto run_dim = third(f0.dim, f1.dim);

      vector3d r;

      r[f0.dim] = bounds_[face0];
      r[f1.dim] = bounds_[face1];
      r[run_dim] = bounds_[run_dim*2];
      auto p0 = r;

      r[run_dim] = bounds_[run_dim*2 + 1];
      auto p1 = r;

      if (
        Face(f0).gridlines_actor_->GetVisibility() &&
        Face(f1).gridlines_actor_->GetVisibility()) {
        // axes_box_.colors_->InsertNextTuple(vector3d{255, 255, 255}); // color per vtk cell, e.g. line
      }
      else {
        axes_box_.colors_->InsertNextTuple(vector3d{255, 255, 255});

        vector2d viewport_coord[2];
      
        vtkNew<vtkCoordinate> coord;
        coord->SetCoordinateSystemToWorld();
        
        coord->SetValue(p0);
        viewport_coord[0] = coord->GetComputedDoubleViewportValue(renderer_);

        // if (viewport_coord[0][0] == 0 && viewport_coord[0][1] == 0) // TODO
        //   return;
        
        coord->SetValue(p1);
        viewport_coord[1] = coord->GetComputedDoubleViewportValue(renderer_);

        if (static constexpr bool flip = cross(vector3i{run_dim}, f0.normal).dot(f1.normal) < 0;
          Face(f0).gridlines_actor_->GetVisibility() ? flip : !flip) {
          std::swap(p0[run_dim], p1[run_dim]);
          std::swap(viewport_coord[0], viewport_coord[1]);
        }
        
        // if (dpl::operation::cross::z(
        //   viewport_coord[1] - viewport_coord[0], 
        //   viewport_coord[2] - viewport_coord[1]) > 0) {
        //
        //   std::swap(p0[run_dim], p1[run_dim]);
        //   std::swap(viewport_coord[0], viewport_coord[1]);
        // }

        // auto viewport_p01 = (viewport_coord[1] - viewport_coord[0]);
        // auto viewport_p01_length = viewport_p01.length();
        auto viewport_axis_dir = (viewport_coord[1] - viewport_coord[0]).normalise(); //viewport_p01/viewport_p01_length;
        auto viewport_tick_dir = vector2d{-viewport_axis_dir.y(), viewport_axis_dir.x()};

        auto alignment = GetAlignment(viewport_tick_dir);

        static constexpr auto tick_length = 5.0;
        static constexpr auto label_to_length = 3.0;

        auto count = count_[run_dim];

        r = p0;

        auto point_idx = axes_ticks_.points_->GetNumberOfPoints();

        int label_skip;

        vector2i label_size;
        
        {
          r[run_dim] = min_[run_dim];
          coord->SetValue(r);
          vector2d t0 = coord->GetComputedDoubleViewportValue(renderer_);

          r[run_dim] += step_[run_dim];
          coord->SetValue(r);
          vector2d t1 = coord->GetComputedDoubleViewportValue(renderer_);

          auto first_tick_step_viewport = t1 - t0;

          r[run_dim] = min_[run_dim] + step_[run_dim]*(count - 2);
          coord->SetValue(r);
          vector2d t2 = coord->GetComputedDoubleViewportValue(renderer_);

          r[run_dim] += step_[run_dim];
          coord->SetValue(r);
          vector2d t3 = coord->GetComputedDoubleViewportValue(renderer_);

          auto last_tick_step_viewport = t3 - t2;

          auto tick_step_viewport = vector2d{
            std::min(std::abs(first_tick_step_viewport.x()), std::abs(last_tick_step_viewport.x())),
            std::min(std::abs(first_tick_step_viewport.y()), std::abs(last_tick_step_viewport.y()))
          };
          
          auto& text_actor = label_actors_[500 - 1]; // last text to test size
          auto* mapper = static_cast<vtkTextMapper*>(text_actor->GetMapper());

          
          
          #ifdef __cpp_lib_format
          {
            double value = (min_[run_dim] + step_[run_dim]*(count - 1))/scale_[run_dim];
            mapper->SetInput(std::vformat(format, std::make_format_args(value)).c_str());
          }
          #else
          {
            double value = (min_[run_dim] + step_[run_dim]*(count - 1))/scale_[run_dim];
            mapper->SetInput(fmt::vformat(format, fmt::make_format_args(value)).c_str());
          }
          // mapper->SetInput((boost::format(format_ptr) %
          //   ((min_[run_dim] + step_[run_dim]*(count - 1))/scale_[run_dim])).str().c_str());
          #endif

          mapper->GetSize(renderer_, label_size);

          label_skip = std::min(
            std::ceil((label_size.x() + 4)/tick_step_viewport.x()),
            std::ceil((label_size.y() + 4)/tick_step_viewport.y()));
        }
        
        for (auto i = 0; i < count; i += label_skip) {
          r[run_dim] = (min_[run_dim] + step_[run_dim]*i)/**scale_[run_dim]*/;
          coord->SetValue(r);
          
          vector2d t0 = coord->GetComputedDoubleViewportValue(renderer_);
          auto t1 = t0 + viewport_tick_dir*tick_length;
          
          axes_ticks_.points_->InsertNextPoint(t0.x(), t0.y(), 0.0);
          axes_ticks_.points_->InsertNextPoint(t1.x(), t1.y(), 0.0);
          
          {
            auto& text_actor = label_actors_[visible_labels_++];
            auto* mapper = static_cast<vtkTextMapper*>(text_actor->GetMapper());

            #ifdef __cpp_lib_format
            {
              double value = r[run_dim]/scale_[run_dim];
              mapper->SetInput(std::vformat(format, std::make_format_args(value)).c_str());
            }
            #else
            {
              double value = r[run_dim]/scale_[run_dim];
              mapper->SetInput(fmt::vformat(format, fmt::make_format_args(value)).c_str());
            }
            // mapper->SetInput((boost::format(format_ptr) % (r[run_dim]/scale_[run_dim])).str().c_str());
            #endif
            
            mapper->GetTextProperty()->SetJustification(alignment.x());
            mapper->GetTextProperty()->SetVerticalJustification(alignment.y());
            
            text_actor->GetPositionCoordinate()->SetValue(
              t1.x() + viewport_tick_dir.x()*label_to_length,
              t1.y() + viewport_tick_dir.y()*label_to_length);
            text_actor->VisibilityOn();
          }
        
          axes_ticks_.lines_->InsertNextCell(2); 
          axes_ticks_.lines_->InsertCellPoint(point_idx++);
          axes_ticks_.lines_->InsertCellPoint(point_idx++);
        }

        {
          auto& text_actor = title_actors_[visible_titles_++];
          auto* mapper = static_cast<vtkTextMapper*>(text_actor->GetMapper());
          mapper->SetInput(titles_[run_dim].c_str());
          mapper->GetTextProperty()->SetJustification(alignment.x());
          mapper->GetTextProperty()->SetVerticalJustification(alignment.y());

          coord->SetValue((p0 + p1)/2);
          vector2d axis_label_pos = coord->GetComputedDoubleViewportValue(renderer_);
        
          vector2d extra_shift{0};

          if (alignment.x() == VTK_TEXT_LEFT)
            extra_shift.x() = label_size.x() /*+ (run_dim == 1 ? 1 : 0)*/;
          else if (alignment.x() == VTK_TEXT_RIGHT)
            extra_shift.x() = -label_size.x() ;

          if (alignment.y() == VTK_TEXT_BOTTOM)
            extra_shift.y() = label_size.y() /*+ (run_dim == 1 ? 1 : 0)*/;
          else if (alignment.y() == VTK_TEXT_TOP)
            extra_shift.y() = -label_size.y();
           
          // auto extra_label_length =
          //   alignment.x() == VTK_TEXT_CENTERED ? /*1.5**/label_size.y() :
          //   alignment.y() == VTK_TEXT_CENTERED ? /*1.5**/label_size.x() :
          //   /*1.2**/std::max(label_size.x(), label_size.y());
          
          
          axis_label_pos =
            axis_label_pos
            + viewport_tick_dir * (tick_length + label_to_length)
            + extra_shift*(
                alignment.x() != VTK_TEXT_CENTERED && alignment.y() != VTK_TEXT_CENTERED ? 1.2 /* Diagonal */ : 
                alignment.x() == VTK_TEXT_CENTERED ? 1.5 /* Top/Bottom */ : 
                1.35 /* Left/Right */
              ); 
        
          text_actor->GetPositionCoordinate()->SetValue(
            axis_label_pos.x(),
            axis_label_pos.y());
          
          text_actor->VisibilityOn();
        }
      }

      
      {
        auto point_idx = axes_box_.points_->GetNumberOfPoints();
                
        axes_box_.points_->InsertNextPoint(p0);
        axes_box_.points_->InsertNextPoint(p1);

        axes_box_.lines_->InsertNextCell(2);
        axes_box_.lines_->InsertCellPoint(point_idx);
        axes_box_.lines_->InsertCellPoint(point_idx + 1);
      }
    }

    template <int face>
    FaceProperties<face>& Face(std::integral_constant<int, face> = {}) {
      return std::get<face>(face_properties_);
    }
    
    vtkNew<vtkActor2D> label_actors_[500];
    int visible_labels_ = 0;

    vector_n<std::string, 3> titles_ = {"X-axis", "Y-axis", "Z-axis"};
    
    vtkNew<vtkActor2D> title_actors_[12];
    int visible_titles_ = 0;

    std::string format_ = ".0f";
  
  public:
    TidyAxes() {
      sfor<6>([&](auto i) {
        Face(i).bounds_ = bounds_.data();
      });

      axes_box_.actor_->UseBoundsOff();
    }

    void SetFormat(const std::string& f) {
      format_ = f;
    }
    
    void Init(
      vtkRenderer* renderer,
      // const StartupSettings& settings,
      const int font_size = 14,
      const char* font_file = nullptr) {
      renderer_ = renderer;

      auto line_width = 1.1f;
        // settings.publishing ? 1.1f : 1.0f;

      vector3d box_ticks_color = {0.6};
      
      {
        axes_box_.poly_data_->SetPoints(axes_box_.points_);
        axes_box_.poly_data_->SetLines(axes_box_.lines_);

        axes_box_.colors_->SetNumberOfComponents(3);
        // axes_box_.poly_data_->GetCellData()->SetScalars(axes_box_.colors_);
       
        axes_box_.mapper_->SetInputData(axes_box_.poly_data_);

        axes_box_.actor_->SetMapper(axes_box_.mapper_);
        axes_box_.actor_->GetProperty()->SetLineWidth(line_width);
        axes_box_.actor_->GetProperty()->SetColor(box_ticks_color);

        
        axes_box_.actor_->UseBoundsOff();
        axes_box_.actor_->PickableOff();
        renderer_->AddActor(axes_box_.actor_);
      }



            
      {

        axes_ticks_.poly_data_->SetPoints(axes_ticks_.points_);
        axes_ticks_.poly_data_->SetLines(axes_ticks_.lines_);

        // axes_ticks_.colors_->SetNumberOfComponents(3);
        // axes_ticks_.poly_data_->GetCellData()->SetScalars(axes_ticks_.colors_);
       
        axes_ticks_.mapper_->SetInputData(axes_ticks_.poly_data_);

        axes_ticks_.actor_->SetMapper(axes_ticks_.mapper_);
        // axes_ticks_.actor_->GetProperty()->width
        axes_ticks_.actor_->GetProperty()->SetLineWidth(line_width);
        axes_ticks_.actor_->GetProperty()->SetColor(box_ticks_color);
        

        axes_ticks_.actor_->UseBoundsOff();
        axes_ticks_.actor_->PickableOff();
        renderer_->AddActor2D(axes_ticks_.actor_);
      }
      



      for (auto& label : label_actors_) {
        vtkNew<vtkTextMapper> text_mapper;
        auto* text_property = text_mapper->GetTextProperty();
        text_property->SetColor(
          vector3d{0.15}
          // settings.publishing ? vector3d{0.15} : vector3d{0}
        );

        if (font_file) {
          text_property->SetFontFamily(VTK_FONT_FILE);
          text_property->SetFontFile(font_file/*, settings.font_path.string().c_str()*/);
        }
        text_property->SetFontSize(font_size/*, settings.visual.font_size*0.85*/);

         // text_property->SetFontFamilyToTimes();
        // text_property->SetFontSize(settings.visual.font_size*0.9/**0.9*1.2*/);
        // text_property->ItalicOn();
        // text_property->BoldOn();
        // text_property->ItalicOn();
        
        label->SetMapper(text_mapper);
        label->GetPositionCoordinate()->SetCoordinateSystemToDisplay();

        renderer_->AddActor2D(label);
        
        label->VisibilityOff();
        label->UseBoundsOff();
        label->PickableOff();

        // label->SetLayerNumber(50);
  
        
      }


      for (auto& label : title_actors_) {
        vtkNew<vtkTextMapper> text_mapper;
        auto* text_property = text_mapper->GetTextProperty();
        text_property->SetColor(
          vector3d{0.15}
          // settings.publishing ? vector3d{0.15} : vector3d{0}
        );
        // text_property->SetFontFamilyToTimes();
        // text_property->SetFontSize(settings.visual.font_size*0.9/**1.2*/);
        // text_property->BoldOn();
        // text_property->ItalicOn();
        // text_property->ShadowOn();
        
        if (font_file) {
          text_property->SetFontFamily(VTK_FONT_FILE);
          text_property->SetFontFile(font_file/*, settings.font_path.string().c_str()*/);
        }
        text_property->SetFontSize(font_size/*settings.visual.font_size*0.85*/);
        
        label->SetMapper(text_mapper);
        label->GetPositionCoordinate()->SetCoordinateSystemToDisplay();

        renderer_->AddActor2D(label);
        
        label->VisibilityOff();
        label->UseBoundsOff();
        label->PickableOff();
      }
      
      
      sfor<6>([&](auto i) {
        Face(i).Init(renderer, vector3d{0.6}, line_width); //settings.publishing ? vector3d{0.6} : vector3d{1}
      });

      renderer_->GetActiveCamera()->AddObserver(vtkCommand::ModifiedEvent, this, &TidyAxes::RefreshAxes);
    }

    void SetScale(const vector3d& scale) {
      scale_ = scale;
    }

    void SetVisibility(bool v) {
      visibility_ = v;
      RefreshAxes();
    }

    bool GetVisibility() const {
      return visibility_;
    }

    

    // void Hide() {
    //   sfor<6>([&](auto i) {
    //     Face(i).Hide();
    //   });
    // }

    void RefreshAxes() {
      sfor<6>([&](auto i) {
        if (visibility_)
          Face(i).RefreshVisibility();
        else
          Face(i).gridlines_actor_->SetVisibility(false);
      });
      
      axes_box_.points_->Initialize();
      axes_box_.lines_->Initialize();
      axes_box_.colors_->Initialize();

      axes_ticks_.points_->Initialize();
      axes_ticks_.lines_->Initialize();

      for (auto i = 0; i < visible_labels_; ++i)
        label_actors_[i]->VisibilityOff();
      visible_labels_ = 0;

      for (auto i = 0; i < visible_titles_; ++i)
        title_actors_[i]->VisibilityOff();
      visible_titles_ = 0;

      sfor<3>([this]<int face0> {
        sfor<face0 + 1, 3>([this]<int face1> {
          RefreshAxis<face0*2, face1*2>();
          RefreshAxis<face0*2, face1*2 + 1>();
          
          RefreshAxis<face0*2 + 1, face1*2>();
          RefreshAxis<face0*2 + 1, face1*2 + 1>();
        });
      });
            
      axes_box_.mapper_->Update();
      axes_ticks_.mapper_->Update();
    }

    void Build(double bounds[6], vector3d min, vector3d step, vector3i count) {
      sfor<6>([&](auto i) {
        bounds_[i] = bounds[i];
      });

      min_ = min;
      step_ = step;
      count_ = count;
      
      sfor<6>([&](auto i) {
        Face(i).BuildGridlines(min, step, count);
      });

      RefreshAxes();
    }

    void Build() {
      Build(renderer_->ComputeVisiblePropBounds());
    }

    /**
     * \param bounds x_min, x_max, y_min, ...
     */
    void Build(double bounds[6]) {
      vector3d min{bounds[0], bounds[2], bounds[4]};

      auto size = vector3d{bounds[1], bounds[3], bounds[5]} - min;
      vector3i count = (size/size.min()*4).min(7.0);
      auto step = size/(count - 1);
      
      min = min/scale_;
      step = step/scale_;
      
      auto adjust_axis = [] (double& min, double& step, int& count) {
        auto value = step;
        auto multipler = std::pow(10.0, std::floor(std::log10(value)));
        
        auto max = (count - 1)*step + min;

        value /= multipler;
        
        value = 
          value < 2.5 ? 1.0 :
          value < 5.0 ? 2.5 :
          5.0;
        
        value *= multipler;

        min = std::ceil(min/value)*value;
        step = value;
        count = (max - min)/step + 1; 
      };

      for (auto i = 0; i < 3; ++i)
        adjust_axis(min[i], step[i], count[i]);
      
      min = min*scale_;
      step = step*scale_;

      constexpr auto margin_coef = 0.075;
      auto flat_factor = std::min({
        bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]})*margin_coef;

      for (auto i = 0; i < 3; ++i) {
        bounds[i*2] -= flat_factor;
        bounds[i*2 + 1] += flat_factor;
      }
      
      Build(bounds, min, step, count);
    }
  };
}





// r[rel.k] = bounds_[rel.k*2];
// coord->SetValue(r);
// face_triangle_display[3] = coord->GetComputedDoubleViewportValue(viewport_);

//
// int on_screen = 0;
// for (int i = 0; i < 4; ++i) {
//   if (0 <= face_triangle_display[i].x() && face_triangle_display[i].x() <= viewport_size.x() &&
//       0 <= face_triangle_display[i].y() && face_triangle_display[i].y() <= viewport_size.y())
//     ++on_screen;
// }

// if (on_screen >= 3)

// else {
//   auto* cmr = viewport_->GetActiveCamera();
//   gridlines_actor_->SetVisibility(
//     (vector3d_map{cmr->GetPosition()} - vector3d_map{cmr->GetFocalPoint()})
//       .normalise().dot(face.normal) > 0);
// }

// if ((face_idx == 0)) {
//   // std::cout << boost::format("%f\n") % 
//   //   vector3d{coord->GetComputedWorldValue(viewport_)};
//   //   face_triangle_display[0];
//
//     auto* cmr = viewport_->GetActiveCamera();
//
//     vector3d_map pos{cmr->GetPosition()};
//   vector3d_map focal{cmr->GetFocalPoint()};
//   
//   std::cout << boost::format("angle = %f; cos = %f;\n")
//   % (std::acos((pos - focal).normalise().dot(face.normal))/3.14159*180.0 - 90)
//   % (pos - focal).normalise().dot(face.normal);
//
//   // gridlines_actor_->SetVisibility(false);
// }