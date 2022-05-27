# could be added to file above
# but not yet required
# lets leave it here

add_executable(test_transform
    test_transform.cc
)

target_include_directories(test_transform
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(test_transform
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(transform test_transform)


add_executable(test_multipoles_planar
  test_multipoles_planar.cc
  multipoles.cc
  # multipoles_extra.cc
  math_comb.cc
)
add_test(multipoles_planar test_multipoles_planar)

target_include_directories(test_multipoles_planar
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(test_multipoles_planar
  # PRIVATE Boost::boost
  #  quadmath
  # gsl
  # gslcblas
  tpsa_lin
    ${Boost_PRG_EXEC_MONITOR_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

target_include_directories(test_multipoles_planar
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
