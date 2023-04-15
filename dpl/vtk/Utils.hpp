/*
 * This file is part of Rapid Reservoir Modelling Software.
 *   | https://rapidreservoir.org/
 *   | https://bitbucket.org/rapidreservoirmodelling/rrm/
 *
 * Copyright (c) 2022
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * RRM is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * RRM is distributed in the hope that it will be useful,                   
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

#include <vtkLookupTable.h>

namespace dpl::vtk
{
  inline void PopulateLut(vtkLookupTable* lut, const std::vector<std::pair<vtkIdType, vector3d>>& entries) {      
    auto prev = entries.begin();      
    auto curr = prev;
    
    while (++curr != entries.end()) {
      double width = curr->first - prev->first;
      for (auto i = prev->first; i < curr->first; ++i) {        
        auto coef = (i - prev->first)/width;          
        auto interp = lerp(prev->second, curr->second, coef);
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
