find_package(Boost 1.52.0 REQUIRED COMPONENTS
  prg_exec_monitor
  unit_test_framework
  program_options
  )

add_subdirectory(gtpsa)
add_subdirectory(tps)

set(thor_scsi_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})


message(STATUS "gtpsa cpp include dir ${gtpsa_cpp_INCLUDE_DIR}")


message(STATUS "thor scsi include dir ${thor_scsi_INCLUDE_DIR}")
add_subdirectory(thor_scsi)
