cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(xpm_project)


find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)
qt_standard_project_setup()

message("Qt6_DIR: " ${Qt6_DIR})
message("Qt6_VERSION: " ${Qt6_VERSION})


find_package(VTK REQUIRED COMPONENTS 
  CommonColor
  CommonCore
  FiltersSources
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

message("VTK_DIR: " ${VTK_DIR})
message("VTK_VERSION: " ${VTK_VERSION})


add_executable(xpm_project
  

  "mywidg.ui"
  "main.cpp"


  #"threedvis.hpp"  
  
  
  
)


target_compile_features(xpm_project PUBLIC cxx_std_20)


target_link_libraries(xpm_project PRIVATE Qt::Core Qt::Gui Qt::Widgets)
if(WIN32)
  add_custom_command(
    TARGET xpm_project POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${Qt6_DIR}/../../$<$<CONFIG:Debug>:debug/>Qt6/plugins/platforms/
      $<TARGET_FILE_DIR:xpm_project>/platforms/
  )

  #add_custom_command(TARGET xpm_project POST_BUILD 
  #  COMMAND Qt6::windeployqt
  #  ARGS $<TARGET_FILE:xpm_project>
  #)
endif()


target_link_libraries(xpm_project PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS xpm_project MODULES ${VTK_LIBRARIES})
