add_executable(test_drift test_drift.cc)
target_include_directories(test_drift
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_drift
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(drift test_drift)

add_executable(test_mpole test_mpole.cc)
target_include_directories(test_mpole
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_mpole
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(mpole test_mpole)

add_executable(test_thin_kick test_thin_kick.cc)
target_include_directories(test_thin_kick
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_thin_kick
  thor_scsi
  thor_scsi_core
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(test_thin_kick test_thin_kick)

add_executable(test_field_kick
  test_field_kick.cc
  )
add_test(field_kick test_field_kick)
target_include_directories(test_field_kick
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_field_kick
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(test_quadrupole
  check_multipole.cc
  test_quadrupole.cc
  )
add_test(quadrupole test_quadrupole)
target_include_directories(test_quadrupole
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_quadrupole
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(test_sextupole
    check_multipole.cc
    test_sextupole.cc
    )
add_test(sextupole test_sextupole)
target_include_directories(test_sextupole
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_sextupole
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(test_cavity test_cavity.cc)
target_include_directories(test_cavity
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_cavity
  thor_scsi
  thor_scsi_core
  tpsa_lin
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(cavity test_cavity)
