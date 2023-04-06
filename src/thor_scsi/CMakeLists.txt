include(../../cmake/gpp_warnings.cmake)

add_subdirectory(core)

add_compile_definitions(
  #   GTSPA_ONLY_OPTIMISED_OPS
  SYNCHROTRON_INTEGRALS
)

if(FLAME_INTERNAL)
    set(flame_INCLUDE_DIR "${flame_SOURCE_DIR}/src/")
    message(STATUS "Using internal flame from subdirectory flame include dir ${flame_INCLUDE_DIR}")
else()
    message(STATUS "using external flame library")
endif()

include_directories(
  ${flame_INCLUDE_DIR}
)


set(thor_scsi_HEADERS
  version.h
  )

# Needs to be fixed ... reorganisation of files
set(thor_scsi_core_HEADERS
  core/exceptions.h
  core/config.h
  core/machine.h
  core/cell_void.h
  core/aperture.h
  core/transform.h
  core/transform_phase_space.h
  core/math_comb.h
  core/field_interpolation.h
  core/multipole_types.h
  core/multipoles.h
  core/elements_basis.h
  core/internals.h
)

set(thor_scsi_element_HEADERS
  elements/constants.h
  elements/drift.h
  elements/elements_enums.h
  elements/utils.h
  elements/radiation_delegate.h
  elements/radiation_delegate_api.h
  elements/field_kick_api.h
  elements/field_kick.h
  elements/element_local_coordinates.h
  elements/mpole.h
  elements/classical_magnet.h
  elements/bending.h
  elements/quadrupole.h
  elements/sextupole.h
  elements/octupole.h
  elements/corrector.h
  elements/element_helpers.h
  elements/cavity.h
  elements/marker.h
  elements/bpm.h
  # Or should that already be a util ...
  elements/standard_observer.h
  elements/standard_aperture.h

  # need g
  )

set(thor_scsi_custom_HEADERS
  custom/aircoil_interpolation.h
  custom/nonlinear_kicker_interpolation.h
  )

set(thor_scsi_std_machine_HEADERS
  std_machine/std_machine.h
  std_machine/accelerator.h
  )

set(thor_scsi_core_FILES
  core/transform.cc
  core/cell_void.cc
  core/elements_basis.cc
  core/machine.cc
  core/math_comb.cc
  core/field_interpolation.cc
  core/multipoles.cc
  core/aperture.cc
  # Only required if GSL's implementation of Horner's rule to be used
  # or a pure taylor series
  # core/multipoles_extra.cc
  elements/utils.cc
  elements/element_helpers.cc
  elements/drift.cc
  elements/cavity.cc
  elements/radiation_delegate.cc
  elements/field_kick.cc
  elements/mpole.cc
  elements/marker.cc
  elements/standard_observer.cc
  elements/standard_aperture.cc
  std_machine/std_machine.cc
  std_machine/accelerator.cc

  custom/aircoil_interpolation.cc
  custom/nonlinear_kicker_interpolation.cc

  # to be removed as soon as tps is removed ...
  ss_vect_tps.cc
  )

add_library(thor_scsi_core SHARED
  ${thor_scsi_core_FILES}
  ${thor_scsi_core_HEADERS}
)
add_library(thor_scsi::core ALIAS thor_scsi_core)

target_include_directories(thor_scsi_core
    PUBLIC
    "$<BUILD_INTERFACE:${thor_scsi_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)


target_link_libraries(thor_scsi_core
  # Todo: need to learn how to instruct cmake that this
  # dependency will be exported properly
  gtpsa-c++
  gtpsa
  tpsa_lin
  flame::core
  # ${flame_CORE_LIBRARY}
  ${ARMADILLO_LIBRARIES}
)

set_target_properties(thor_scsi_core
  PROPERTIES
    SOVERSION 0
    LIBRARY_OUTPUT_NAME thor_scsi_core_gtpsa
)


set(thor_scsi_files
  version.cc
)

add_library(thor_scsi SHARED
  ${thor_scsi_files}
  ${thor_scsi_HEADERS}
)

message(STATUS "gtpsa_cpp_INCLUDEDIR:"  ${gtpsa_cpp_INCLUDE_DIR})

target_include_directories(thor_scsi
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<BUILD_INTERFACE:${gtpsa_cpp_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set_target_properties(thor_scsi
  PROPERTIES
    SOVERSION 0
    LIBRARY_OUTPUT_NAME thor_scsi_gtpsa
)

add_executable(test_transform_phase_space
  core/test_transform_phase_space.cc
)
target_link_libraries(test_transform_phase_space
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

add_test(transform_phase_space test_transform_phase_space)

## add_executable(a_kick  elements/a_kick.cc)
## target_include_directories(a_kick
##     PUBLIC
##     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
##     "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
## )
## target_link_libraries(a_kick
##   thor_scsi
##   thor_scsi_core
##   tpsa_lin
##   # PRIVATE Boost::boost
##   #  quadmath
##     ${Boost_PRG_EXEC_MONITOR_LIBRARY}
##     ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
## )

#
# Non linear part needs still to be ported to new interface
## set(tpsa_nlin_FILES
##     tpsa_for_pm.cc
##     # tpsa_for.cc
## )
##
## add_library(tpsa_nlin SHARED
##   ${tpsa_common_FILES}
##   ${tpsa_nlin_FILES}
##   ${tpsa_HEADERS}
## )



## set(thor_scsi_CORE_HEADERS
##   thor/core/util.h
##   flame/core/base.h
##   flame/core/config.h
##   )
##
##
##
## add_library(thor_scsi SHARED
##   ${thor_scsi_CORE_files}
##   ${thor_scsi_CORE_HEADERS}
## )


add_subdirectory(elements)
add_subdirectory(std_machine)
add_subdirectory(examples)
add_subdirectory(custom)

# ---- install helpers --

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(THOR_SCSI_CMAKE_LIB_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/thor_scsi")

message(STATUS "CMAKE install Dirs  lib: '${CMAKE_INSTALL_LIBDIR}' bin '${CMAKE_INSTALL_BINDIR}' include '${CMAKE_INSTALL_INCLUDEDIR}'")


# --- thor_scsi thor install support  ----

install(FILES ${thor_scsi_core_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
install(FILES ${thor_scsi_element_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/elements/)
install(FILES ${thor_scsi_std_machine_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/std_machine/)
install(FILES ${thor_scsi_custom_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/custom/)
install(TARGETS thor_scsi_core
        EXPORT thor_scsi_coreTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
# generate and install export file
install(EXPORT thor_scsi_coreTargets
        FILE coreTargets.cmake
        NAMESPACE thor_scsi::
	DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)

# --- end thor_scsi thor install support -

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/thor_scsiConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfig.cmake"
    INSTALL_DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfigVersion.cmake
    VERSION "${THOR_SCSI_VERSION}"
    COMPATIBILITY SameMajorVersion
)


# --- thor scsi install support  ---------
install(FILES ${thor_scsi_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/)
install(FILES ${thor_scsi_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/)
install(TARGETS thor_scsi
        EXPORT thor_scsiTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT thor_scsiTargets
    FILE thor_scsiTargets.cmake
    DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
    NAMESPACE thor_scsi::
)

# --- end thor scsi install support  -----

# --- tpsa install support  -----

install(FILES
         "${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfig.cmake"
         "${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfigVersion.cmake"
	 DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)
