
if(CMAKE_COMPILER_IS_GNUCXX)
  message(STATUS "Enabling gxx debug preprocessor flags")
  # add_definitions(-Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers)
  add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g -ggdb -pedantic -Wall)
endif()

# yes ... depends on flame ...
find_package(flame  REQUIRED COMPONENTS core)
get_target_property(flame_CORE_LIBRARY flame::flame_core IMPORTED_LOCATION_RELEASE)
message(STATUS "Flame include dir ${flame_INCLUDE_DIR}")
get_filename_component(thor_scsi_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}//../inc/" ABSOLUTE)

set(tpsa_HEADERS
    ${thor_scsi_INCLUDE_DIR}/tps/ss_vect.h
    ${thor_scsi_INCLUDE_DIR}/tps/mmf.h
    ${thor_scsi_INCLUDE_DIR}/tps/tps_type.h
    ${thor_scsi_INCLUDE_DIR}/tps/math_pass.h
    ${thor_scsi_INCLUDE_DIR}/tps/enums.h
    ${thor_scsi_INCLUDE_DIR}/tps/externs.h
    ${thor_scsi_INCLUDE_DIR}/tps/forward_decl.h
    ${thor_scsi_INCLUDE_DIR}/tps/tpsa_for_pm.h
    ${thor_scsi_INCLUDE_DIR}/tps/tpsa_for.h
    ${thor_scsi_INCLUDE_DIR}/tps/tpsa_lin_pm.h
#    ${thor_scsi_INCLUDE_DIR}/tps/field.h
    ${thor_scsi_INCLUDE_DIR}/tps/config.h
    ${thor_scsi_INCLUDE_DIR}/tps/exceptions.h
    ${thor_scsi_INCLUDE_DIR}/tps/tpsa.h
    ${thor_scsi_INCLUDE_DIR}/tps/tps.h
    ${thor_scsi_INCLUDE_DIR}/tps/utils.h
    ${thor_scsi_INCLUDE_DIR}/tps/ss_vect_utils.h
    ${thor_scsi_INCLUDE_DIR}/tps/tpsa_lin.h
)


set(tpsa_common_FILES
    tpsa_utils.cc
    field.cc
)

set(tpsa_lin_FILES
    tpsa_lin.cc
    tpsa_lin_pm.cc
)


add_library(tpsa_lin SHARED
  ${tpsa_common_FILES}
  ${tpsa_lin_FILES}
  ${tpsa_HEADERS}
)


get_filename_component(thor_scsi_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../inc/" ABSOLUTE)
message(STATUS "TPS(A) include directory ${thor_scsi_TPS_INCLUDE_DIR}")

include_directories(
  ${thor_scsi_INCLUDE_DIR}
  ${flame_INCLUDE_DIR}
)

add_library(thor_scsi::tpsa_lin ALIAS tpsa_lin)

target_include_directories(tpsa_lin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../inc/>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)


set(thor_scsi_HEADERS
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/version.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/exceptions.h
  )

# Needs to be fixed ... reorganisation of files
set(thor_scsi_CORE_HEADERS
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/transform.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/lattice_filter.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/cells.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/constants.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements_basis.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/defines.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/core.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements_enums.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/elements_filter.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/element_helpers.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/lattice.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/config.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/cpp_version.h
)
set(thor_scsi_CORE_ELEMENT_HEADERS
  elements/drift.h
)

set(thor_scsi_FILES
  core/transform.cc
  elements/drift.cc
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

add_library(thor SHARED
  ${thor_scsi_FILES}
  ${thor_scsi_HEADERS}
  ${thor_scsi_CORE_HEADERS}
)
add_library(thor_scsi::thor ALIAS thor)

target_include_directories(thor
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)


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

# ---- install helpers --
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(THOR_SCSI_CMAKE_LIB_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/thor_scsi")
# --- tpsa install support  --------------
install(FILES ${tpsa_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tps/)
install(TARGETS tpsa_lin
        EXPORT thor_scsitpsa_linTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
# generate and install export file
install(EXPORT thor_scsitpsa_linTargets
        FILE thor_scsitpsa_linTargets.cmake
        NAMESPACE thor_scsi::
	DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)
# --- end tpsa install support  ----------
# --- thor_scsi install support  ---------
install(FILES ${thor_scsi_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/)
install(FILES ${thor_scsi_CORE_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
install(FILES ${thor_scsi_CORE_ELEMENT_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
install(TARGETS thor
        EXPORT thor_scsithorTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
# generate and install export file
install(EXPORT thor_scsithorTargets
        FILE thor_scsithorTargets.cmake
        NAMESPACE thor_scsi::
	DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)

# --- end thor_scsi install support  -----

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
# --- tpsa install support  -----

install(FILES
         "${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfig.cmake"
         "${CMAKE_CURRENT_BINARY_DIR}/thor_scsiConfigVersion.cmake"
	 DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)
