cmake_minimum_required(VERSION 3.0.0)

project("thor_scsi" LANGUAGES CXX Fortran)

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
  set(CMAKE_BUILD_TYPE RELEASE) # default to debug build
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
set(mad_ng_SRC_DIR ${gtpsa_DIR}/mad-ng/src/)

include_directories(
  ${mad_ng_SRC_DIR}/libgtpsa/
  ${mad_ng_SRC_DIR}
  )

# Fortran implementation part
add_subdirectory(TPSA)
add_subdirectory(src)
# python wrapper
add_subdirectory(python)
