cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0148)
  # allow legacy FindPythonInterp/FindPythonLibs modules
  cmake_policy(SET CMP0148 OLD)
  set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
endif()
if(POLICY CMP0167)
  # keep classic FindBoost behavior to avoid config-package issues
  cmake_policy(SET CMP0167 OLD)
  set(CMAKE_POLICY_DEFAULT_CMP0167 OLD)
endif()

project("thor_scsi" LANGUAGES CXX Fortran)

set(Boost_NO_WARN_NEW ON)
set(Boost_NO_BOOST_CMAKE OFF)
set(Boost_NO_WARN_NEW ON CACHE BOOL "Silence Boost new-version warnings" FORCE)
set(Boost_NO_BOOST_CMAKE OFF CACHE BOOL "Prefer FindBoost over config packages" FORCE)
set(Boost_ADDITIONAL_VERSIONS "1.89.0" "1.89")
set(Boost_ADDITIONAL_VERSIONS "1.89.0;1.89" CACHE STRING "Known Boost versions" FORCE)

# OpenMP (AppleClang needs libomp from Homebrew)
set(OpenMP_ROOT "/opt/homebrew/opt/libomp" CACHE PATH "OpenMP root path")
set(OPENMP_ROOT "/opt/homebrew/opt/libomp" CACHE PATH "OpenMP root path (upper)")
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp")
if(APPLE AND CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
  set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "" FORCE)
  set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
  set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib" CACHE FILEPATH "" FORCE)
  include_directories("/opt/homebrew/opt/libomp/include")
endif()

# MacOS - set RPATH entries.
set(CMAKE_BUILD_RPATH "${CMAKE_SOURCE_DIR}/local/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/local/lib")

# Use absolute vs. relative paths.
# For macOS - RPATH & DYLD_LIBRARY.
# set(CMAKE_INSTALL_NAME_DIR "${CMAKE_SOURCE_DIR}/local/lib")
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH)
# set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/local/lib")
# New in version 3.9.
# set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)

ENABLE_LANGUAGE(Fortran)

set(THOR_SCSI_VERSION 0.0.1)

# option(PYTHON_MODULE_INSTALL "Install python module" OFF)
option(PYTHON_MODULE_DEVELOPMENT_INSTALL "Install python binary module to source directory" ON)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(-g -ggdb -pedantic -Wall)
  # add_definitions(
  #  -Wall
  #  -Wextra
  #  -Wno-unused-parameter
  #  -Wno-missing-field-initializers
  #  )
  #  message(STATUS "Enabling gxx debug preprocessor flags")
  # flame relies on boost ... have to search if ubuntu has a version which supports
  # GNU extra flags ...
  #  add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG) # default to debug build
  # set(CMAKE_BUILD_TYPE RELEASE) # default to debug build
  message(STATUS "Default to DEBUG build")
endif()

find_package(Armadillo REQUIRED)

find_package(PythonInterp 3)
# find_package(Java REQUIRED)
# message(STATUS "Java include dirs" ${JAVA_INCLUDE_PATH})
# find_package(JNI REQUIRED)
# message(STATUS "Jini include dirs" ${JINI_INCLUDE_DIRS})
# include(CMakeAddFortranSubdirectory)

enable_testing()

set(FLAME_INTERNAL True)
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


# Need to write proper .cmake
set(gtpsa_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/gtpsa/)
set(gtpsa_cpp_INCLUDE_DIR ${gtpsa_DIR}/c++)
# only required by test transfomr?
set(gtpsa_mad_ng_INCLUDE_DIR ${gtpsa_DIR}/mad-ng/src)

# Fortran implementation part
#add_subdirectory(TPSA)
add_subdirectory(src)
# python wrapper
add_subdirectory(python)
