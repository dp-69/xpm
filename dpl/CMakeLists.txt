#
# This file is part of Dmytro Petrovskyy Library (DPL).
# 
# Copyright (c) 2023
#   | Dmytro Petrovskyy, PhD
#   | dmytro.petrovsky@gmail.com
#   | https://www.linkedin.com/in/dmytro-petrovskyy/
# 
# DPL is free software: you can redistribute it and/or modify              
# it under the terms of the GNU General Public License as published by     
# the Free Software Foundation, either version 3 of the License, or        
# (at your option) any later version.                                      
#                                                                         
# DPL is distributed in the hope that it will be useful,                   
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
# GNU General Public License for more details.                             
#                                                                         
# You should have received a copy of the GNU General Public License        
# along with RRM. If not, see <http://www.gnu.org/licenses/>.
# 
#

add_library(dpl INTERFACE)

target_include_directories(dpl INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")

option(DPL_HYPRE_MPI "<help_text>")

if(DPL_HYPRE_MPI)
  target_compile_definitions(dpl INTERFACE DPL_HYPRE_BOOST_SHARED_MEMORY)
endif()