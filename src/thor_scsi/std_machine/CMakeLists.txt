add_executable(test_accelerator test_accelerator.cc)

target_include_directories(test_accelerator
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_accelerator
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(accelerator test_accelerator)

add_executable(test_accelerator_config test_accelerator_config.cc)

target_include_directories(test_accelerator_config
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_accelerator_config
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(accelerator test_accelerator_config)
