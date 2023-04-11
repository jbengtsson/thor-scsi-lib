add_executable(test_aircoil_interpolation test_aircoil_interpolation.cc)
target_include_directories(test_aircoil_interpolation
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_aircoil_interpolation
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa-c++
  gtpsa
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(aircoil_interpolation test_aircoil_interpolation)

add_executable(test_nonlinear_kicker_interpolation test_nonlinear_kicker_interpolation.cc)
target_include_directories(test_nonlinear_kicker_interpolation
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_nonlinear_kicker_interpolation
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa-c++
  gtpsa
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(nonlinear_kicker_interpolation test_nonlinear_kicker_interpolation)

add_executable(test_nonlinear_kicker test_nonlinear_kicker.cc)
target_include_directories(test_nonlinear_kicker
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_nonlinear_kicker
  thor_scsi
  thor_scsi_core
  tpsa_lin
  gtpsa-c++
  gtpsa
  # PRIVATE Boost::boost
  #  quadmath
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(nonlinear_kicker test_nonlinear_kicker)
