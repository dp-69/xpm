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

#include <vtkLookupTable.h>

namespace dpl::vtk
{
  struct CameraSettings
  {
    bool parallel = false;
    double parallel_scale = 1.0;
    double view_angle;
    vector3d position;
    vector3d focal_point;
    vector3d view_up;

    CameraSettings() = default;

    explicit CameraSettings(vtkCamera* cam) {
      parallel = cam->GetParallelProjection();
      parallel_scale = cam->GetParallelScale();
      cam->GetPosition(position);
      cam->GetViewUp(view_up);
      cam->GetFocalPoint(focal_point);
      view_angle = cam->GetViewAngle();
    }

    void operator>>(vtkCamera* cam) {
      cam->SetParallelProjection(parallel);
      cam->SetParallelScale(parallel_scale);
      cam->SetPosition(position);
      cam->SetViewUp(view_up);
      cam->SetFocalPoint(focal_point);
      cam->SetViewAngle(view_angle);
    }
  };

  inline void to_json(nlohmann::json& j, const CameraSettings& p) {
    j = nlohmann::json{
      {"parallel",       p.parallel},
      {"parallel_scale", p.parallel_scale},
      {"position",       p.position},
      {"focal_point",    p.focal_point},
      {"view_up",        p.view_up},
      {"view_angle",     p.view_angle}
    };
  }
  
  inline void from_json(const nlohmann::json& j, CameraSettings& p) {
    p.position       = j.at("position");
    p.focal_point    = j.at("focal_point");
    p.view_up        = j.at("view_up");
    p.parallel       = j.at("parallel");
    p.parallel_scale = j.at("parallel_scale");
    p.view_angle     = j.at("view_angle");
  }


  inline void PopulateLut(vtkLookupTable* lut, const std::vector<std::pair<vtkIdType, vector3d>>& entries) {      
    auto prev = entries.begin();      
    auto curr = prev;
    
    while (++curr != entries.end()) {
      double width = curr->first - prev->first;               // NOLINT(clang-diagnostic-implicit-int-float-conversion, bugprone-narrowing-conversions, cppcoreguidelines-narrowing-conversions)
      for (auto i = prev->first; i < curr->first; ++i) {          
        auto coef = (i - prev->first)/width;                  // NOLINT(clang-diagnostic-implicit-int-float-conversion, bugprone-narrowing-conversions)
        auto interp = dpl::lerp(prev->second, curr->second, coef);
        lut->SetTableValue(i, interp.x(), interp.y(), interp.z());
      }
        
      prev = curr;
    }

    lut->SetTableValue(prev->first, prev->second.x(), prev->second.y(), prev->second.z());           
  }

  inline void PopulateLutRedWhiteBlue(vtkLookupTable* lut, const int entries = 1000) {
    lut->SetNumberOfTableValues(entries);
    
    PopulateLut(lut, {
      {0, {0.231373, 0.298039, 0.752941}},
      {3*entries/10, {0.548188, 0.581521, 0.808972}}, //
      {entries/2, {0.865003, 0.865003, 0.865003}},
      {7*entries/10, {0.785443, 0.440345, 0.507012}}, //
      {entries - 1, {0.705882, 0.0156863, 0.14902}}
    });
  }
}
