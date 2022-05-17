find_package(pybind11 2.9 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)


# detect virtualenv and set Pip args accordingly
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()



## # install PyPI Python package using pip
## execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${_pip_args} ${_pypkg_name})
## # install Python package locally in development mode
## execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${_pip_args} -e .)

cmake_policy(SET CMP0063 NEW)
add_library(lib MODULE
  src/thor_scsi.cc
  src/arma.cc
  src/tps.cc
  src/config_type.cc
  src/elements.cc
  src/accelerator.cc
  )

message(STATUS "python wrapper flame include dir ${flame_INCLUDE_DIR}")
target_include_directories(lib
  PUBLIC
  ${ARMADILLO_INCLUDE_DIRS}
  ${flame_INCLUDE_DIR}
  # Hack for pybind11 on Pierre's laptop
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../flame/src/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(lib
  pybind11::module
  pybind11::lto
  pybind11::windows_extras
  thor_scsi
  thor_scsi_core
  tpsa_lin
  ${flame_CORE_LIBRARY}
  ${ARMADILLO_LIBRARIES}
  #
  #
  )

pybind11_extension(lib)
pybind11_strip(lib)
set_target_properties(lib PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")

add_library(pyflame MODULE
  src/flame.cc
  )

target_include_directories(pyflame
  PUBLIC
  ${flame_INCLUDE_DIR}
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../flame/src/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(pyflame
  pybind11::module
  pybind11::lto
  pybind11::windows_extras
  thor_scsi
  thor_scsi_core
  tpsa_lin
  ${flame_CORE_LIBRARY}
  #
  #
  )
pybind11_extension(pyflame)
pybind11_strip(pyflame)
set_target_properties(pyflame
  PROPERTIES
  CXX_VISIBILITY_PRESET "hidden"
  CUDA_VISIBILITY_PRESET "hidden"
  )


# For development ... let's copy the built shared objects to thor_scsi
# For final installation: use setup.py

##add_custom_target(copy-thor-scsi-lib ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib.cpython-38-x86_64-linux-gnu.so)
##add_custom_command(
##                   OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/thor_scsi/lib.cpython-38-x86_64-linux-gnu.so
##                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/lib.cpython-38-x86_64-linux-gnu.so ${CMAKE_CURRENT_SOURCE_DIR}/thor_scsi
##                   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib.cpython-38-x86_64-linux-gnu.so)
