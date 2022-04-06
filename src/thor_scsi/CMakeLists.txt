include(../../cmake/gpp_warnings.cmake)

# yes ... depends on flame ...
find_package(flame  REQUIRED COMPONENTS core)
get_target_property(flame_CORE_LIBRARY flame::flame_core IMPORTED_LOCATION_RELEASE)



message(STATUS "Flame core library ${flame_CORE_LIBRARY}")
message(STATUS "Flame include dir ${flame_INCLUDE_DIR}")


# get_filename_component(thor_scsi_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}//../inc/" ABSOLUTE)
# get_filename_component(thor_scsi_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../inc/" ABSOLUTE)

add_subdirectory(core)

include_directories(
#  ${thor_scsi_INCLUDE_DIR}
  ${flame_INCLUDE_DIR}
)


set(thor_scsi_HEADERS
  version.h
  )

# Needs to be fixed ... reorganisation of files
set(thor_scsi_core_HEADERS
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/exceptions.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/lattice_filter.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/constants.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/defines.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/core.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements_enums.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements_filter.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/lattice.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/config.h
  # ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/cpp_version.h
  core/exceptions.h
  core/config.h
  core/machine.h
  core/cell_void.h
  core/cells.h
  core/transform.h
  core/transform_phase_space.h
  core/math_comb.h
  core/field_interpolation.h
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
  elements/element_helpers.h
  elements/cavity.h
  elements/marker.h
  elements/bpm.h
  # Or should that already be a util ...
  elements/standard_observer.h
  )

set(thor_scsi_std_machine_HEADERS
  std_machine/std_machine.h
  std_machine/accelerator.h
)
set(thor_scsi_core_FILES
  core/transform.cc
  core/cell_void.cc
  core/machine.cc
  core/math_comb.cc
  core/field_interpolation.cc
  core/multipoles.cc
  # Only required if GSL's implementation of Horner's rule to be used
  # or a pure taylor series
  # core/multipoles_extra.cc
  elements/utils.cc
  elements/element_helpers.cc
  elements/drift.cc
  elements/radiation_delegate.cc
  elements/field_kick.cc
  elements/mpole.cc
  elements/marker.cc
  elements/cavity.cc
  elements/standard_observer.cc
  std_machine/std_machine.cc
  std_machine/accelerator.cc
  # thor_lib.cc
  # t2cell.cc
  # t2lat.cc
  # radia2tracy.cc
  # legacy/io.cc
  # ety.cc
  # eigenv.cc
  # prtmfile.cc
  # rdmfile.cc
  # lsoc.cc
  # set_errors.cc
)

add_library(thor_scsi_core SHARED
  ${thor_scsi_core_FILES}
  ${thor_scsi_core_HEADERS}
)
add_library(thor_scsi::core ALIAS thor_scsi_core)

target_include_directories(thor_scsi_core
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)


target_link_libraries(thor_scsi_core
  # Todo: need to learn how to instruct cmake that this
  # dependency will be exported properly
  tpsa_lin
  ${flame_CORE_LIBRARY}
  ${ARMADILLO_LIBRARIES}
)

set_target_properties(thor_scsi_core
  PROPERTIES
    SOVERSION 0
)


set(thor_scsi_files
  version.cc
)

add_library(thor_scsi SHARED
  ${thor_scsi_files}
  ${thor_scsi_HEADERS}
)

target_include_directories(thor_scsi
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)


add_executable(test_transform_phase_space
  core/test_transform_phase_space.cc
)
target_link_libraries(test_transform_phase_space
  thor_scsi
  thor_scsi_core
  tpsa_lin
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

# ---- install helpers --

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(THOR_SCSI_CMAKE_LIB_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/thor_scsi")

message(STATUS "CMAKE install Dirs  lib: '${CMAKE_INSTALL_LIBDIR}' bin '${CMAKE_INSTALL_BINDIR}' include '${CMAKE_INSTALL_INCLUDEDIR}'")


# --- thor_scsi thor install support  ----

install(FILES ${thor_scsi_core_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
install(FILES ${thor_scsi_element_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/elements/)
install(FILES ${thor_scsi_std_machine_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/std_machine/)
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
