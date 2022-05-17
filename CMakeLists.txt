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

find_package(Armadillo REQUIRED)

find_package(PythonInterp 3)

# include(CMakeAddFortranSubdirectory)

enable_testing()

set(FLAME_INTERNAL TRUE)
if(FLAME_INTERNAL)
    add_subdirectory(flame)
else()
    # yes ... depends on flame ...
    find_package(flame  REQUIRED COMPONENTS core)

    # Flame used by t
    get_target_property(flame_CORE_LIBRARY flame::flame_core IMPORTED_LOCATION_RELEASE)
    message(STATUS "Flame core library ${flame_CORE_LIBRARY}")
    message(STATUS "Flame include dir ${flame_INCLUDE_DIR}")
endif()

# Fortran implementation part
add_subdirectory(TPSA)
# C++ implementation part
add_subdirectory(src)
# add_subdirectory(thor)
add_subdirectory(python)
