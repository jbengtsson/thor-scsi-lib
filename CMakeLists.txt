cmake_minimum_required(VERSION 3.0.0)

project("thor_scsi" LANGUAGES CXX Fortran)
# project("thor_scsi")

ENABLE_LANGUAGE(Fortran)

set(THOR_SCSI_VERSION 0.0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(-g -ggdb -pedantic -Wall)
  add_definitions(-Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers)
  #  message(STATUS "Enabling gxx debug preprocessor flags")
  # flame relies on boost ... have to search if ubuntu has a version which supports
  # GNU extra flags ...
  #  add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE) # default to debug build
  message(STATUS "Default to DEBUG build")
endif()


find_package(PythonInterp 3)

# include(CMakeAddFortranSubdirectory)

add_subdirectory(TPSA)
add_subdirectory(thor)
