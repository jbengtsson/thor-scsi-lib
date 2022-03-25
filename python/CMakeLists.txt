find_package(pybind11 2.9 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

#pybind11_add_module(flame
#  src/flame.cc
# )


cmake_policy(SET CMP0063 NEW)
add_library(lib MODULE
  src/thor_scsi.cc
  src/tps.cc
  src/config_type.cc
  src/elements.cc
  src/machine.cc
  )

target_include_directories(lib
  PUBLIC
  # Hack for pybind11 on Pierre's laptop
    "$<BUILD_INTERFACE:/home/mfp/.local/lib/python3.8/site-packages/pybind11/include/>"
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
  #
  #
  )

pybind11_extension(lib)
pybind11_strip(lib)
set_target_properties(lib PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")

add_library(flame MODULE
  src/flame.cc
  )

target_include_directories(flame
  PUBLIC
  # Hack for pybind11 on Pierre's laptop
    "$<BUILD_INTERFACE:/home/mfp/.local/lib/python3.8/site-packages/pybind11/include/>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../flame/src/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(flame
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
pybind11_extension(flame)
pybind11_strip(flame)
set_target_properties(flame PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")
