find_package(Boost 1.52.0 REQUIRED COMPONENTS
  prg_exec_monitor
  unit_test_framework
  program_options
  )

add_subdirectory(tps)
add_subdirectory(thor_scsi)
