# -*-  mode: CMAKE; cmake-tab-width: 4 -*-
# $Id: CMakeLists.txt 384 2009-08-24 21:04:09Z doudoubai $

INCLUDE_DIRECTORIES(${soap_SOURCE_DIR}/glps)
ADD_LIBRARY(glps parser scanner)

# Create target for the parser
IF(MSVC)
  #message(STATUS "Using MS VS ")
  SET(CMAKE_CXX_FLAGS "/MT")
  SET(CMAKE_C_FLAGS "/MT")
    # hdf5
  #INCLUDE_DIRECTORIES(${DEVEL_ROOT}/hdf5/include)
ENDIF(MSVC)


FIND_PROGRAM(LEX_EXE flex)
IF(LEX_EXE STREQUAL "LEX_EXE-NOTFOUND")
    MESSAGE(FATAL_ERROR "ERROR: flex not found!")
ENDIF(LEX_EXE STREQUAL "LEX_EXE-NOTFOUND")

FIND_PROGRAM(YACC_EXE bison )
IF(YACC_EXE STREQUAL "YACC_EXE-NOTFOUND")
  MESSAGE(FATAL_ERROR "ERROR: yacc not found!")
ENDIF(YACC_EXE STREQUAL "YACC_EXE-NOTFOUND")

MESSAGE(STATUS "Found lex : ${LEX_EXE}")
MESSAGE(STATUS "Found yacc: ${YACC_EXE}")

# Create custom command for flex/lex (note the outputs)
ADD_CUSTOM_COMMAND(
  SOURCE ${soap_SOURCE_DIR}/glps/scanner.l
  COMMAND ${LEX_EXE} 
  ARGS -i -o${soap_BINARY_DIR}/glps/scanner.c ${soap_SOURCE_DIR}/glps/scanner.l
  TARGET dak
  OUTPUTS ${soap_BINARY_DIR}/glps/scanner.c)

# Create custom command for bison/yacc (note the DEPENDS)
ADD_CUSTOM_COMMAND(
  SOURCE ${soap_SOURCE_DIR}/glps/parser.y
  COMMAND ${YACC_EXE}
  ARGS -t -v -d -y ${soap_SOURCE_DIR}/glps/parser.y -o ${soap_BINARY_DIR}/glps/parser.c
  TARGET dak
  DEPENDS ${soap_BINARY_DIR}/glps/scanner.c
  OUTPUTS ${soap_BINARY_DIR}/glps/parser.c)
#endif(WIN32)
# Add parser.c to the list of sources
#SET(dak_SRCS ${dak_SRCS} ${dak_BINARY_DIR}/src/parser.c beamline.cc latio.cc)

# Since parser.c does not exists yet when cmake is run, mark
# it as generated
SET_SOURCE_FILES_PROPERTIES(${soap_BINARY_DIR}/glps/parser.c GENERATED)
SET_SOURCE_FILES_PROPERTIES(${soap_BINARY_DIR}/glps/parser.h GENERATED)
SET_SOURCE_FILES_PROPERTIES(${soap_BINARY_DIR}/glps/scanner.c GENERATED)

#IF(MSVC)
#    SET_SOURCE_FILES_PROPERTIES(${soap_BINARY_DIR}/src/parser.c 
#	COMPILE_FLAGS "/Tp")
#ENDIF(MSVC)

# Include binary directory to include lexer.c in parser.c
#INCLUDE_DIRECTORIES(${soap_BINARY_DIR}/src)

#
#ADD_EXECUTABLE(latxml latxml)

