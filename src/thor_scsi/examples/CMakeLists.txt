add_executable(lattice_pass lattice_pass.cc)

target_include_directories(lattice_pass
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(lattice_pass
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa
  Boost::program_options
  #  quadmath
)

add_executable(build_accelerator build_accelerator.cc)

target_include_directories(build_accelerator
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(build_accelerator
  thor_scsi
  thor_scsi_core
  tpsa_lin
  #  quadmath
)
