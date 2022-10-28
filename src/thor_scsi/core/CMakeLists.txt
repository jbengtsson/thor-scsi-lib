# could be added to file above
# but not yet required
# lets leave it here

add_executable(test_transform
    test_transform.cc
)

message(STATUS "thor_scsi::core gtpsa cpp include dir ${gtpsa_cpp_INCLUDE_DIR}")
target_include_directories(test_transform
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
    "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(test_transform
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(transform test_transform)


add_executable(test_two_dimensional_multipoles
  test_two_dimensional_multipoles.cc
  multipoles.cc
  # multipoles_extra.cc
  math_comb.cc
)
add_test(two_dimensional_multipoles test_two_dimensional_multipoles)

target_include_directories(test_two_dimensional_multipoles
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
    "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_two_dimensional_multipoles
  # PRIVATE Boost::boost
  #  quadmath
  # gsl
  # gslcblas
  tpsa_lin
  gtpsa-c++
  gtpsa
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
