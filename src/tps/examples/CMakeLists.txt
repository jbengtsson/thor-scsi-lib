add_executable(a_mul_tps a_mul.cc)
target_include_directories(a_mul_tps
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
# "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"

target_link_libraries(a_mul_tps
  tpsa_lin
  )

add_executable(a_ss_vect a_ss_vect.cc)
target_include_directories(a_ss_vect
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
# "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"

target_link_libraries(a_ss_vect
  tpsa_lin
  )
