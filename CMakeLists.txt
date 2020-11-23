cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(Liggghts_VBP_VTKReader)

find_package(VTK COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkIOLegacy QUIET)
find_package(Python2 COMPONENTS Development)
if (NOT VTK_FOUND)
  message("Skipping Liggghts_VBP_VTKReader: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

set(SOURCE_FILES
        src/Liggghts_VBP_VTKReader.cxx
        src/sphericalParticle.cpp src/superquadricParticle.cpp src/particle.cpp src/bedHeight.h)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(Liggghts_VBP_VTKReader MACOSX_BUNDLE ${SOURCE_FILES})
  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES} )
#  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES} Python2::Python)
#  target_include_directories(Liggghts_VBP_VTKReader PRIVATE ${Python2_INCLUDE_DIRS})
  if (WIN32)
    target_include_directories(Liggghts_VBP_VTKReader PRIVATE include)
  endif (WIN32)
else ()
  # include all components
  add_executable(Liggghts_VBP_VTKReader MACOSX_BUNDLE ${SOURCE_FILES})
  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES})
#  target_link_libraries(Liggghts_VBP_VTKReader PRIVATE ${VTK_LIBRARIES} Python2::Python)
#  target_include_directories(Liggghts_VBP_VTKReader PRIVATE ${Python2_INCLUDE_DIRS})
  if (WIN32)
    target_include_directories(Liggghts_VBP_VTKReader PRIVATE include)
  endif (WIN32)
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS Liggghts_VBP_VTKReader
    MODULES ${VTK_LIBRARIES}
    )
endif ()
