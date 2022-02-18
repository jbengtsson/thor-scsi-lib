cmake_minimum_required(VERSION 3.0.0)

project("thor_scsi" LANGUAGES CXX Fortran)
# project("thor_scsi")

ENABLE_LANGUAGE(Fortran)

set(THOR_SCSI_VERSION 0.0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(PythonInterp 3)

# include(CMakeAddFortranSubdirectory)

add_subdirectory(TPSA)
add_subdirectory(thor)
