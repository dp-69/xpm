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

#include <QVariant>
#include <QWidget>

namespace dpl::qt::property_editor
{
  template <typename Scalar> struct ScalarConvert {};

  template<> struct ScalarConvert<bool> {
    static bool FromString(const QString& str, bool& value) {
      value = QVariant{str}.toBool();
      return true;
    }

    static auto ToQString(bool value) { return QString::number(value); }
    static auto ToString(bool value) { return std::to_string(value); } 
  };
  
  template<> struct ScalarConvert<double> {
    static bool FromString(const QString& str, double& value) {
      bool valid;
      auto str_copy = str;
      value = str_copy.replace(",", "").toDouble(&valid);
      return valid;
    }

    static auto ToQString(double value) { return QString::number(value); }
    static auto ToString(double value) { return std::to_string(value); } 
  };

  template<> struct ScalarConvert<int> {
    static bool FromString(const QString& str, int& value) {
      bool valid;
      value = str.toInt(&valid);
      return valid;
    }
    
    static auto ToQString(int value) { return QString::number(value); }
    static auto ToString(int value) { return std::to_string(value); } 
  };

  template<> struct ScalarConvert<QString> {
    static auto FromString(const QString& str, QString& value) {
      value = str;
      return true;
    }
    
    static auto ToQString(const QString& value) { return value; }
    static auto ToString(const QString& value) { return value.toStdString(); }
  };

  template<> struct ScalarConvert<std::string> {
    static auto FromString(const QString& str, std::string& value) {
      value = str.toStdString();
      return true;
    }

    static auto ToQString(const std::string& value) { return QString::fromStdString(value); }
    static auto ToString(const std::string& value) { return value; }
  };
}