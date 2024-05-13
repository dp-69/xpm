#
# This file is part of Dmytro Petrovskyy Library (dpl).
# 
# Copyright (c) 2023
#   | Dmytro Petrovskyy, PhD
#   | dmytro.petrovsky@gmail.com
#   | https://www.linkedin.com/in/dmytro-petrovskyy/
# 
# dpl is free software: you can redistribute it and/or modify              
# it under the terms of the GNU General Public License as published by     
# the Free Software Foundation, either version 3 of the License, or        
# (at your option) any later version.                                      
#                                                                         
# dpl is distributed in the hope that it will be useful,                   
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
# GNU General Public License for more details.                             
#                                                                         
# You should have received a copy of the GNU General Public License        
# along with dpl. If not, see <http://www.gnu.org/licenses/>.
# 
#

add_library(dpl INTERFACE)

target_include_directories(dpl INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")