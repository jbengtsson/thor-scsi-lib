
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    message(STATUS "Compiling with GNU Fortran ${CMAKE_Fortran_COMPILER_VERSION}")
    # Inconsistencies in size of integers detected by compiler on MSYS platform
    # I guess due to g++ version 11: now moving to default
    add_compile_options(-g -O2 -Wall -fbounds-check
      -fdefault-integer-8 -finteger-4-integer-8
      -w
      -mcmodel=medium
      )
    add_compile_options(-fpermissive)
    #    add_compile_options(-fcheck-array-temporaries)
    add_compile_options(-ff2c)
    if(MSYS)
        # Inconsistencies in size of integers detected by compiler on MSYS platform
        # I guess due to g++ version 11
        # perhaps also due to platform
        # should be fixed in source code
        # add_compile_options(-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large)
        # add_compile_options(-g
        #                    -O2
    	# 		    -Wall
    	#		    -fbounds-check
    	#		    -fdefault-integer-8
	#		    -finteger-4-integer-8
        #		    -w
    	#	             -mcmodel=medium
    	#		    )
        # add_compile_options(-fpermissive)
      endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
    add_compile_options(-g -O2 -Wall -fbounds-check
      -fdefault-integer-8 -finteger-4-integer-8
      -w
      -mcmodel=medium)
    add_compile_options(-fpermissive)
    add_compile_options(-ff2c)
    message(STATUS "Compiling with flang")
else()
  message(STATUS "Not compiling with GNU Fortran nor with flang but ${CMAKE_Fortran_COMPILER_ID}")
endif()

set(tpsa_lie_files ./LieLib.f)
set(tpsa_tpsa_files ./TPSALib.f)

# set(CMAKE_Fortran_CREATE_SHARED_LIBRARY ar)

if(MSYS)
    # As gfortran-11 will not compile it
    # should  be a check for the compiler version
    message(STATUS "Not compiling LieLIb")
else()
    message(INFO "TPSA lie files ${tpsa_lie_files} tpsa files ${tpsa_tpsa_files}")
    add_library(LieLib STATIC  ${tpsa_lie_files})
    set_target_properties(LieLib PROPERTIES LINKER_LANGUAGE Fortran)
endif()

add_library(TPSALib STATIC  ${tpsa_tpsa_files})
set_target_properties(TPSALib PROPERTIES LINKER_LANGUAGE Fortran)
