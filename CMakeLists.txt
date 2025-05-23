#
# This file is part of Extensive Pore Modelling (xpm).
#   | https://github.com/dp-69/xpm
#
# Copyright (c) 2024
#   | Dmytro Petrovskyy, PhD
#   | dmytro.petrovsky@gmail.com
#   | https://www.linkedin.com/in/dmytro-petrovskyy/
#
# xpm is free software: you can redistribute it and/or modify              
# it under the terms of the GNU General Public License as published by     
# the Free Software Foundation, either version 3 of the License, or        
# (at your option) any later version.                                      
#                                                                         
# xpm is distributed in the hope that it will be useful,                   
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
# GNU General Public License for more details.                             
#                                                                         
# You should have received a copy of the GNU General Public License        
# along with xpm. If not, see <http://www.gnu.org/licenses/>.
#
# 

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(xpm_project)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_AUTOMOC ON)

if(MSVC)  
  string(REPLACE "/W3" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST:NO")

  #set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /MANIFEST:NO")
  #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /MANIFEST:NO")
endif()

find_package(Qt6 COMPONENTS Core Gui Widgets OpenGLWidgets Charts)

if(Qt6_FOUND)
  qt_standard_project_setup()
  message("Qt6_DIR: " ${Qt6_DIR})
  message("Qt6_VERSION: " ${Qt6_VERSION})
endif()


find_package(VTK REQUIRED)

message("VTK_DIR: " ${VTK_DIR})
message("VTK_VERSION: " ${VTK_VERSION})


add_executable(xpm # WIN32
  "src/widget.hpp"  
  "src/main.cpp"
  
  
  "libs/dpl/qt/property_editor/Editors.hpp"
  "libs/dpl/qt/property_editor/PropertyEntities.hpp"
  "libs/dpl/qt/property_editor/QPropertyTreeView.hpp"
)

target_include_directories(xpm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_compile_features(xpm PUBLIC cxx_std_20)

target_compile_definitions(xpm PRIVATE XPM_DEBUG_OUTPUT)


add_subdirectory(libs/dpl)
target_link_libraries(xpm PRIVATE dpl)

add_subdirectory(libs/pnextract)

if(MSVC)  
  target_compile_options(pnextract PRIVATE -w)
  target_compile_options(xpm PRIVATE -W4)

  #if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  #  target_compile_options(xpm PRIVATE -bigobj)
  #endif()

  target_compile_options(xpm PRIVATE -bigobj)
else()
  #add_compile_options("-w")
endif()

set_target_properties(pnextract
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pnextract")

if(Qt6_FOUND)
  target_link_libraries(xpm PRIVATE Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGLWidgets Qt6::Charts)

  if(WIN32)
    add_custom_command(
      TARGET xpm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${Qt6_DIR}/../../$<$<CONFIG:Debug>:debug/>Qt6/plugins/platforms/
        $<TARGET_FILE_DIR:xpm>/platforms/
    )

    add_custom_command(
      TARGET xpm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${Qt6_DIR}/../../$<$<CONFIG:Debug>:debug/>Qt6/plugins/styles/
        $<TARGET_FILE_DIR:xpm>/styles/
    )
  endif()
endif()

target_link_libraries(xpm PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS xpm MODULES VTK::RenderingFreeType) # Needed in Ubuntu

find_package(Boost COMPONENTS iostreams graph REQUIRED)
target_link_libraries(xpm PRIVATE Boost::iostreams Boost::graph)

find_package(argh CONFIG REQUIRED)
target_link_libraries(xpm PRIVATE argh)

find_package(fmt CONFIG REQUIRED)
target_link_libraries(xpm PRIVATE fmt::fmt)

find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(xpm PRIVATE nlohmann_json::nlohmann_json)


find_package(MPI)
if(MPI_FOUND)
  target_link_libraries(xpm PRIVATE MPI::MPI_CXX)
endif()


if(TRUE)
  set(HYPRE_ENABLE_HYPRE_BLAS OFF)
  set(HYPRE_ENABLE_HYPRE_LAPACK OFF)
  set(HYPRE_ENABLE_MIXEDINT OFF)
  set(HYPRE_ENABLE_BIGINT OFF)
  set(HYPRE_ENABLE_SINGLE OFF)
  set(HYPRE_PRINT_ERRORS OFF)
  add_subdirectory("libs/hypre/src")
  target_link_libraries(xpm PRIVATE HYPRE)
else()
  find_package(HYPRE REQUIRED)
  target_link_libraries(xpm PRIVATE HYPRE::HYPRE)
endif()


if(NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.json)
  file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  file(COPY_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/files/config.json"
    "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.json")
endif()