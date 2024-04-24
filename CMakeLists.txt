# Minimum required Version
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
# Generates a 'compile_commands.json' file containing the exact compiler calls for all translation units of the project in machine-readable form.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Since GNU Compiler is always available, look first for mpiifort which is usually linked with intel
#find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} mpiifort mpif90 PATHS ENV PATH NO_DEFAULT_PATH)  # If ever using MPI, use this line instead
find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} ifort gfortran PATHS ENV PATH NO_DEFAULT_PATH)

# Project name and languages (C or CXX must be activated to findLAPACK from MKL)
project(PROJECT Fortran C)

# Defining options
option(DEBUG "Enable debugging" OFF)
option(PREP "Enable Instrumentation" OFF)
#find_package(MPI) # If ever using MPI, use this line

# Source files
set(SOURCES     source/mod_global.f90
                source/bogoliubov_transf.f90
                source/spinwavesNonCol.f90
                )

# Getting the COMPILER when not given
if(NOT COMPILER)
  set(COMPILER ${CMAKE_Fortran_COMPILER_ID})
  message("Setting COMPILER=${COMPILER}")
endif()
# Changing $COMPILER variable to lowercase for comparison later
string(TOLOWER ${COMPILER} COMPILER)

# Adding debugging suffix when present
if(DEBUG)
  set(DEBUGSUFFIX _debug)
  message("DEBUG on")
else()
  set(DEBUGSUFFIX)
endif()

# Adding compiler suffix:
# INTEL
if(COMPILER MATCHES intel)
  set(CMAKE_Fortran_COMPILER mpiifort)
  set(COMPILERSUFFIX _intel)
  message("Using Intel Compiler")
endif()
# GNU
if(COMPILER MATCHES gfortran OR COMPILER MATCHES gcc OR COMPILER MATCHES gfort OR COMPILER MATCHES gnu)
  set(COMPILER gnu)
  set(CMAKE_Fortran_COMPILER mpif90)
  set(COMPILERSUFFIX _gnu)
  message("Using GNU Compiler")
endif()

# Adding other suffixes when present
if(SUFFIX)
  set(SUFFIX _${SUFFIX})
else()
  set(SUFFIX)
endif()

######################################## EXECUTABLE ########################################

# Defining executable name
#set(PROJECT_EXE main${DEBUGSUFFIX}${COMPILERSUFFIX}${INSTRUMENT}${COVSUFFIX}${SUFFIX}.exe)
set(PROJECT_EXE main.exe)

# Add an executable to the project using the specified source files
add_executable(${PROJECT_EXE} ${SOURCES})
# Set (main) language
set_property(TARGET ${PROJECT_EXE} PROPERTY LINKER_LANGUAGE Fortran)
# Set location for binary
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

######################################## LIBRARIES ########################################
# Add libraries
find_package(BLAS)

find_package(LAPACK REQUIRED)
find_package(OpenMP)
#target_link_libraries(${PROJECT_EXE} MPI::MPI_Fortran LAPACK::LAPACK OpenMP::OpenMP_Fortran) 
target_link_libraries(${PROJECT_EXE} LAPACK::LAPACK OpenMP::OpenMP_Fortran)
#target_link_libraries(${PROJECT_EXE} OpenMP::OpenMP_Fortran)


# Getting GIT version
execute_process(COMMAND git describe --abbrev=8 --dirty --always OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

# Setting DEBUG definition
if(DEBUG)
  target_compile_definitions(${PROJECT_EXE} PRIVATE DEBUG)
endif(DEBUG)

######################################## COMPILATION FLAGS ########################################

if(COMPILER MATCHES intel) # INTEL:
  if(NOT DEBUG)
    target_compile_options(${PROJECT_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -qopt-zmm-usage=high
                                                -align array64byte
                                                -march=core-avx2
                                                -O3
                                                )
  else(NOT DEBUG)
    target_compile_options(${PROJECT_EXE} PRIVATE -C
                                                -CB
                                                -check uninit
                                                -debug all
                                                -warn all
                                                -e18
                                                -ftrapuv
                                                -traceback
                                                -fpp
                                                -fpe0
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -O0
                                                )
  endif(NOT DEBUG)
elseif(COMPILER MATCHES gnu) # GFORTRAN:
  if(NOT DEBUG)
    target_compile_options(${PROJECT_EXE} PRIVATE -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -fallow-argument-mismatch
                                                -ffree-form
                                                -ffree-line-length-0
                                                -g
                                                -O3
                                                -march=native
                                                )

  else(NOT DEBUG)
    target_compile_options(${PROJECT_EXE} PRIVATE -C
                                                -Wall
                                                -Wextra
                                                -Wconversion
                                                -Wline-truncation
                                                -Wcharacter-truncation
                                                -Wsurprising
                                                -Waliasing
                                                -Wunused-parameter
                                                -Wno-maybe-uninitialized
                                                -fwhole-file
                                                -fimplicit-none
                                                -fallow-argument-mismatch
                                                -fbacktrace
                                                -ffree-form
                                                -ffree-line-length-0
                                                -fcheck=all
                                                #-ffpe-trap=zero,overflow,underflow
                                                -finit-real=nan
                                                -ftrapv
                                                -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -O0
                                                )
  endif(NOT DEBUG)
endif(COMPILER MATCHES intel)

# unset(DEBUG CACHE)
# unset(DEBUGSUFFIX CACHE)
# unset(PREP CACHE)
# unset(SUFFIX CACHE)
# unset(COMPILER CACHE)
