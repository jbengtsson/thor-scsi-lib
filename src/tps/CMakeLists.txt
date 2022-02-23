# Add g++ warnings for development
include(../../cmake/gpp_warnings.cmake)

set(tpsa_HEADERS
    ss_vect.h
    mmf.h
    tps_type.h
    math_pass.h
    enums.h
    externs.h
    forward_decl.h
    tpsa_for_pm.h
    tpsa_for.h
    tpsa_lin_pm.h
#   field.h
    config.h
    exceptions.h
    tpsa.h
    tps.h
    utils.h
    ss_vect_utils.h
    tpsa_lin.h
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
add_library(thor_scsi::tpsa_lin ALIAS tpsa_lin)

target_include_directories(tpsa_lin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)

set_target_properties(tpsa_lin
  PROPERTIES
    SOVERSION 0
)

target_link_libraries(tpsa_lin
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --- tpsa install support  --------------
install(FILES ${tpsa_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tps/)
install(TARGETS tpsa_lin
        EXPORT tpsa_linTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
