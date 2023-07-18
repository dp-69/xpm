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

inline dpl::vector3d hsb_to_rgb(double hue, double sat, double bright) {
  if (hue > 360 || hue < 0 || sat > 100 || sat < 0 || bright > 100 || bright < 0) {
    // cout<<"The givem HSV values are not in valid range"<<endl;
    return {};
  }
  auto s = sat/100;
  auto v = bright/100;
  auto C = s*v;
  auto X = C*(1 - std::abs(std::fmod(hue/60.0, 2) - 1));
  auto m = v - C;

  double r, g, b;
  if (hue >= 0 && hue < 60) {
    r = C;
    g = X;
    b = 0;
  }
  else if (hue >= 60 && hue < 120) {
    r = X;
    g = C;
    b = 0;
  }
  else if (hue >= 120 && hue < 180) {
    r = 0;
    g = C;
    b = X;
  }
  else if (hue >= 180 && hue < 240) {
    r = 0;
    g = X;
    b = C;
  }
  else if (hue >= 240 && hue < 300) {
    r = X;
    g = 0;
    b = C;
  }
  else {
    r = C;
    g = 0;
    b = X;
  }

  return {r, g, b};
  
  // int R = (r + m) * 255;
  // int G = (g + m) * 255;
  // int B = (b + m) * 255;
}
