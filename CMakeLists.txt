cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(Liggghts_VBP_VTKReader)

find_package(VTK COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkIOLegacy QUIET)

if (NOT VTK_FOUND)
  message("Skipping Liggghts_VBP_VTKReader: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(Liggghts_VBP_VTKReader MACOSX_BUNDLE Liggghts_VBP_VTKReader.cxx)
  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES} )
  target_include_directories(Liggghts_VBP_VTKReader PRIVATE)
else ()
  # include all components
  add_executable(Liggghts_VBP_VTKReader MACOSX_BUNDLE Liggghts_VBP_VTKReader.cxx)
  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES})
  target_include_directories(Liggghts_VBP_VTKReader PRIVATE )
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS Liggghts_VBP_VTKReader
    MODULES ${VTK_LIBRARIES}
    )
endif ()
