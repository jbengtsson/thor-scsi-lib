find_package(pybind11 2.9 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

cmake_policy(SET CMP0063 NEW)
add_library(lib MODULE
  src/thor_scsi.cc
  src/tps.cc
  src/config_type.cc
  src/aperture.cc
  src/interpolation.cc
  src/radiation.cc
  src/observer.cc
  src/elements.cc
  src/accelerator.cc
  )

message(STATUS "python wrapper flame include dir ${flame_INCLUDE_DIR}")
target_include_directories(lib
  PUBLIC
  ${ARMADILLO_INCLUDE_DIRS}
  ${flame_INCLUDE_DIR}
  "$<BUILD_INTERFACE:${thor_scsi_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${flame_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(lib
  pybind11::module
  pybind11::lto
  pybind11::windows_extras
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa
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
    "$<BUILD_INTERFACE:${thor_scsi_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${flame_INCLUDE_DIR}>"
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
if(PYTHON_MODULE_DEVELOPMENT_INSTALL)
  message(STATUS "copying python dynamic built libraries to source dir ${CMAKE_CURRENT_SOURCE_DIR}")
  install(TARGETS pyflame DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/thor_scsi)
  install(TARGETS lib DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/thor_scsi)
endif()

# detect virtualenv and set Pip args accordingly
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()


if(PYTHON_MODULE_INSTALL)
    # not proper code (nor tested)
    # execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install ${_pip_args} ${_pypkg_name})
endif()
if(PYTHON_MODULE_DEVELOPMENT_INSTALL)
    # not proper code (nor tested)
    # install Python package locally in development mode
    #  execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install ${_pip_args} -e )
endif()
