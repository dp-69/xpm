#
# This file is part of Rapid Reservoir Modelling Software.
#   | https://rapidreservoir.org/
#   | https://bitbucket.org/rapidreservoirmodelling/rrm/
#
# Copyright (c) 2022
#   | Dmytro Petrovskyy, PhD
#   | dmytro.petrovsky@gmail.com
#   | https://www.linkedin.com/in/dmytro-petrovskyy/
#
# RRM is free software: you can redistribute it and/or modify              
# it under the terms of the GNU General Public License as published by     
# the Free Software Foundation, either version 3 of the License, or        
# (at your option) any later version.                                      
#                                                                         
# RRM is distributed in the hope that it will be useful,                   
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