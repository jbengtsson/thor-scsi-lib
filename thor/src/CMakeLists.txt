
# Hope that spills not into fortran code
if(CMAKE_COMPILER_IS_GNUCXX)

  add_definitions(
  	-Wcast-align
	-Wcast-qual
	-Wctor-dtor-privacy
	-Wdisabled-optimization
	-Wformat=2
	-Winit-self
	-Wlogical-op
	-Wnoexcept
	# -Wold-style-cast
	-Woverloaded-virtual
	-Wredundant-decls
	# -Wshadow
	-Wsign-promo
	-Wstrict-null-sentinel
	-Wstrict-overflow=5
	-Wswitch-default
	-Wundef
	-Wno-unused
   	-Wmissing-declarations
   	#-Weffc++
	#-Wpadding
)
  # has to be removed at a later stage
  # guess a move ctor and use in code will provide
  # some speed up
  # add_definitions(-Werror)

  add_definitions(
  -Wno-deprecated-copy
    -Wno-error=deprecated-copy
  )
  add_definitions(
  -Wno-Werror=all
   -Wmissing-declarations
  -Wno-error=unused-result
  -Wno-error=sign-compare
  -Wno-error=vla -Wno-error=old-style-cast
  -Wno-error=shadow
  -Wno-error=strict-overflow
  -Wno-error=sign-promo
  -Wno-error=format-nonliteral
  -Wno-error=switch-default)
endif()

# yes ... depends on flame ...
find_package(flame  REQUIRED COMPONENTS core)
get_target_property(flame_CORE_LIBRARY flame::flame_core IMPORTED_LOCATION_RELEASE)
message(STATUS "Flame core library ${flame_CORE_LIBRARY}")
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

set_target_properties(tpsa_lin
  PROPERTIES
    SOVERSION 0
)


target_link_libraries(tpsa_lin
)



set(thor_scsi_HEADERS
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/version.h
  )

# Needs to be fixed ... reorganisation of files
set(thor_scsi_core_HEADERS
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/exceptions.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/transform.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/lattice_filter.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/cell_void.h
  ${thor_scsi_INCLUDE_DIR}/thor_scsi/core/machine.h
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
set(thor_scsi_core_ELEMENT_HEADERS
  elements/drift.h
)

set(thor_scsi_core_FILES
  core/transform.cc
  core/cell_void.cc
  core/machine.cc
  core/element_helpers.cc
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

add_library(thor_scsi_core SHARED
  ${thor_scsi_core_FILES}
  ${thor_scsi_core_HEADERS}
)
add_library(thor_scsi::core ALIAS thor_scsi_core)

target_include_directories(thor_scsi_core
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)


target_link_libraries(thor_scsi_core
# Todo: need to learn how to include it ...
#  tpsa_lin
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

message(STATUS "CMAKE install Dirs  lib: '${CMAKE_INSTALL_LIBDIR}' bin '${CMAKE_INSTALL_BINDIR}' include '${CMAKE_INSTALL_INCLUDEDIR}'")


# --- tpsa install support  --------------
install(FILES ${tpsa_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tps/)
install(TARGETS tpsa_lin
        EXPORT tpsa_linTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
# generate and install export file
install(EXPORT tpsa_linTargets
        FILE tpsa_linTargets.cmake
        NAMESPACE thor_scsi::
	DESTINATION ${THOR_SCSI_CMAKE_LIB_DIR}
)
# --- end tpsa install support  ----------

# --- thor_scsi thor install support  ----

install(FILES ${thor_scsi_core_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
install(FILES ${thor_scsi_core_ELEMENT_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thor_scsi/core/)
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
