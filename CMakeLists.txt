cmake_minimum_required(VERSION 3.0.1)
project(Droplets_2D)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

#~ find_package(VTK REQUIRED)
#~ include(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/solver)

#ADD_SUBDIRECTORY(src/cost_function_library)
#ADD_SUBDIRECTORY(src/local_optimization_library)

ADD_EXECUTABLE(Droplets_2D src/main.cpp)

ADD_DEPENDENCIES(Droplets_2D options ini_parser utils solver)

TARGET_LINK_LIBRARIES(Droplets_2D options ini_parser utils solver dl m)
