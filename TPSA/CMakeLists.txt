
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  message(STATUS "Compiling with GNU Fortran ${CMAKE_Fortran_COMPILER_VERSION}")
  add_definitions(-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium -w -O2)
else()
  message(STATUS "Not compiling with GNU Fortran")
endif()

set(tpsa_lie_files ./LieLib.f)
set(tpsa_tpsa_files ./TPSALib.f)

# set(CMAKE_Fortran_CREATE_SHARED_LIBRARY ar)

message(STATUS "TPSA lie files ${tpsa_lie_files} tpsa files ${tpsa_tpsa_files}")
add_library(LieLib STATIC  ${tpsa_lie_files})

set_target_properties(LieLib PROPERTIES LINKER_LANGUAGE Fortran)

add_library(TPSALib STATIC  ${tpsa_tpsa_files})
set_target_properties(TPSALib PROPERTIES LINKER_LANGUAGE Fortran)
