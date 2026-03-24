# (C) Copyright 2025 Meteorolgisk Institutt
#

# ECSABER bundle
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.6 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)
project( ecsaber-bundle VERSION 0.0.1 LANGUAGES C CXX Fortran )
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

# Options
option( ENABLE_MPI "Compile with MPI" ON )
option( ENABLE_OMP "Compile with OpenMP" ON )
option( ENABLE_CUDA "Compile with CUDA" OFF )
option( ENABLE_MKL "Compile with MKL" OFF )
option( ENABLE_TRANSI "Compile TransI C-interface to trans" ON )
option( ENABLE_ETRANS "Include Limited-Area-Model Transforms" ON )
option( ENABLE_JEDI_MODEL_DATA "Obtain saber/vader test data from jedi-model-data repository (vs tarball)" ON )
option( ENABLE_QUENCH "Compile quench" ON )

# Activate ECSABER
set( ECSABER ON CACHE BOOL "Compile for ECSABER" )

# Initialize bundle
ecbuild_bundle_initialize()

# Use external jedi-cmake or build in bundle
if( DEFINED ENV{jedi_cmake_ROOT} )
    include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
    list( APPEND CMAKE_MODULE_PATH $ENV{jedi_cmake_ROOT}/share/jedicmake/Modules )
    message( INFO "CMAKE_MODULE_PATH" )
else()
    ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda/jedi-cmake.git" BRANCH develop UPDATE RECURSIVE )
    include( jedicmake/cmake/Functions/git_functions.cmake )
    list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/jedicmake/CMakeModules/Modules )
endif()

# ECMWF tools
ecbuild_bundle( PROJECT eckit    GIT "https://github.com/ecmwf/eckit.git" TAG 1.26.3 )
ecbuild_bundle( PROJECT fckit    GIT "https://github.com/ecmwf/fckit.git" TAG 0.11.2 )
ecbuild_bundle( PROJECT fiat     GIT "https://github.com/hirlam/fiat.git " BRANCH fiat-CY49T2h_oops )
ecbuild_bundle( PROJECT ectrans  GIT "https://github.com/accord-nwp/ectrans.git" BRANCH for-CY49T1 )
ecbuild_bundle( PROJECT atlas    GIT "https://github.com/benjaminmenetrier/atlas.git" BRANCH atlas-CY49T2h_oops )

# OOPS
#ecbuild_bundle( PROJECT oops GIT "https://github.com/benjaminmenetrier/oops.git" BRANCH feature/CY49T2h_ecsaber_to_CY50_4dvar UPDATE )
ecbuild_bundle( PROJECT oops )

# jedi-data
ecbuild_bundle( PROJECT jedi-model-data GIT "https://github.com/jcsda-internal/jedi-model-data.git" BRANCH develop UPDATE )

# ECSABER
ecbuild_bundle( PROJECT ecsaber SOURCE ..)
set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_BINARY_DIR}/ecsaber/oops_patch )
set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_BINARY_DIR}/ecsaber/saber )
set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_BINARY_DIR}/ecsaber/vader )

ecbuild_bundle_finalize()
