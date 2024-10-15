# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# QUENCHXX_TEST
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( quenchxx_test VERSION 0.0.1 LANGUAGES C CXX Fortran )

option(OPENMP "Build quenchxx_test with OpenMP support" ON)
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})
set( CMAKE_DIRECTORY_LABELS "quenchxx_test" )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

################################################################################
# Dependencies
################################################################################

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS C Fortran )

# atlas
find_package( atlas 0.35.0 REQUIRED )

# oops
find_package( oops 0.0.0 REQUIRED )

# oops_patch
find_package( oops_patch 0.0.1 REQUIRED )

# saber
find_package( saber 0.0.1 REQUIRED )

#OpenMP
if(OPENMP)
  find_package( OpenMP COMPONENTS Fortran )
endif()

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )

################################################################################
# Sources
################################################################################

add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME quenchxx_test )

# print the summary of the configuration
ecbuild_print_summary()
