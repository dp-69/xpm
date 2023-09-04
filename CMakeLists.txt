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
else()
  find_package(Qt5 COMPONENTS Core Gui Widgets Charts REQUIRED) #Charts
  message("Qt5_DIR: " ${Qt5_DIR})
  message("Qt5_VERSION: " ${Qt5_VERSION})
endif()


find_package(VTK REQUIRED)

#find_package(VTK REQUIRED COMPONENTS 
#  CommonColor
#  CommonCore
#  #FiltersSources
#  #InteractionStyle
#  #RenderingContextOpenGL2
#  #RenderingCore
#  #RenderingFreeType
#  #RenderingGL2PSOpenGL2
#  #RenderingOpenGL2
#  GUISupportQt
#)

message("VTK_DIR: " ${VTK_DIR})
message("VTK_VERSION: " ${VTK_VERSION})


add_executable(xpm # WIN32
  "src/xpm_widget.hpp"  
  "src/main.cpp"
  "src/functions.cpp"
  
  "dpl/qt/property_editor/Editors.hpp"
  "dpl/qt/property_editor/PropertyEntities.hpp"
  "dpl/qt/property_editor/QPropertyTreeView.hpp"
)

target_include_directories(xpm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_compile_features(xpm PUBLIC cxx_std_20)

target_compile_definitions(xpm PRIVATE XPM_DEBUG_OUTPUT)


add_subdirectory(dpl)
target_link_libraries(xpm PRIVATE dpl)

add_subdirectory(pnextract)

if(MSVC)  
  target_compile_options(pnextract PRIVATE -w)
  target_compile_options(xpm PRIVATE -W4)
else()
  #add_compile_options("-w")
endif()

set_target_properties(pnextract PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pnextract")

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
else()
  target_link_libraries(xpm PUBLIC Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Charts)
endif()



if(VTK_VERSION VERSION_LESS "8.90.0")
  message("VTK_USE_FILE: " ${VTK_USE_FILE})
  include(${VTK_USE_FILE})

  target_link_libraries(xpm PRIVATE ${VTK_LIBRARIES}
    vtkCommonColor
    vtkCommonCore
    vtkFiltersModeling
    vtkFiltersSources
    vtkChartsCore
    vtkImagingCore
    vtkIOImage
    vtkInteractionStyle
    vtkInteractionWidgets
    vtkRenderingAnnotation    
    vtkRenderingContextOpenGL2
    vtkRenderingCore
    vtkRenderingFreeType
    vtkRenderingGL2PSOpenGL2
    vtkRenderingOpenGL2
    vtkRenderingVolumeOpenGL2
    vtkGUISupportQt)
else()
  target_link_libraries(xpm PRIVATE ${VTK_LIBRARIES})
  vtk_module_autoinit(TARGETS xpm MODULES ${VTK_LIBRARIES})
endif()



find_package(Boost COMPONENTS iostreams graph REQUIRED)
target_link_libraries(xpm PRIVATE Boost::iostreams Boost::graph)


find_package(fmt CONFIG REQUIRED)
target_link_libraries(xpm PRIVATE fmt::fmt)

find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(xpm PRIVATE nlohmann_json::nlohmann_json)


find_package(MPI)
if(MPI_FOUND)
  target_link_libraries(xpm PRIVATE MPI::MPI_CXX)
endif()

find_package(HYPRE REQUIRED)
target_link_libraries(xpm PRIVATE HYPRE::HYPRE)

if(NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.json)
  file(COPY_FILE
    ${CMAKE_CURRENT_SOURCE_DIR}/config/config.json
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.json)
endif()