# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( ecsaber VERSION 0.0.1 LANGUAGES CXX Fortran )

# Update flag
set( ECSABER_UPDATE OFF )
if( DEFINED ENV{ECSABER_UPDATE} )
    set( ECSABER_UPDATE $ENV{ECSABER_UPDATE} )
endif()

# Process ECSABER tools
add_subdirectory( tools )

# Check JCSDA base repos
if( ECSABER_UPDATE )
    if( NOT EXISTS "${CMAKE_SOURCE_DIR}/ecsaber/oops-jedi")
        message( FATAL_ERROR "Missing oops-jedi repository, do: git clone OOPS_URL oops-jedi" )
    endif()
    if( NOT EXISTS "${CMAKE_SOURCE_DIR}/ecsaber/saber-jedi")
        message( FATAL_ERROR "Missing saber-jedi repository, do: git clone SABER_URL saber-jedi" )
    endif()
    if( NOT EXISTS "${CMAKE_SOURCE_DIR}/ecsaber/vader-jedi")
        message( FATAL_ERROR "Missing vader-jedi repository, do: git clone VADER_URL vader-jedi" )
    endif()
endif()

# Update functions
function( patch_init )
    if( ECSABER_UPDATE )
        # Remove update script
        file( REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )

        # Create update script
        file( WRITE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh "#!/usr/bin/env bash\n" )
    endif()
endfunction()

function( check_diff inputFile outputFile )
    if( ECSABER_UPDATE )
        # Run check_diff script
        execute_process( COMMAND bash
                                 ${CMAKE_BINARY_DIR}/bin/ecsaber_check_diff.sh
                                 ${inputFile}
                                 ${outputFile}
                                 "ECSABER"
                                 ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )
    endif()
endfunction()

function( check_diff_recursively outputDir )
    if( ECSABER_UPDATE )
        # Output directory in ECSABER
        set( ECSABER_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${outputDir})

        # Input directory in JEDI
        string( REPLACE "ecsaber/${PROJECT_NAME}" "ecsaber/${PROJECT_NAME}-jedi" JEDI_SOURCE_DIR ${ECSABER_SRC_DIR} )

        if( IS_DIRECTORY ${ECSABER_SRC_DIR} )
            # Run check_diff_recursively script
            execute_process( COMMAND bash
                                     ${CMAKE_BINARY_DIR}/bin/ecsaber_check_diff_recursively.sh
                                     ${CMAKE_BINARY_DIR}/bin/ecsaber_check_diff.sh
                                     ${JEDI_SOURCE_DIR}
                                     ${ECSABER_SRC_DIR}
                                     "ECSABER"
                                     ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )
        else()
            # Run check_diff script
            execute_process( COMMAND bash
                                     ${CMAKE_BINARY_DIR}/bin/ecsaber_check_diff.sh
                                     ${JEDI_SOURCE_DIR}
                                     ${ECSABER_SRC_DIR}
                                     "ECSABER"
                                     ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )
        endif()
    endif()
endfunction()

function( patch_apply )
    if( ECSABER_UPDATE )
        # Run update script
        execute_process( COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh
                         COMMAND_ERROR_IS_FATAL ANY )

        # Remove update script
        file( REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )
    endif()
endfunction()

# Download SABER test files
message( STATUS "Checking SABER test files" )
set( saber_test_files
CovStats.nc
dirac_gsi_geos_global.nml
dirac_gsi_geos_global_opt_1.nml
dirac_gsi_geos_global_opt_2.nml
dirac_gsi_geos_global_opt_3.nml
dirac_gsi_gfs_global.nml
FPstats.nc
gauss_state.nc
gsi-coeffs-global-l127x64y32.nc4
gsi-coeffs-gmao-global-l72x72y46.nc4
gsi-coeffs-gmao-global-l72x72y46_opt1.nc4
gsi-coeffs-gmao-global-l72x72y46_opt2.nc4
gsi-coeffs-gmao-global-l72x72y46_opt3.nc4
Lv.nc
MIO_coefficients.nc
MUstats.nc
Prho_bar_Mean.nc
ptheta_bar_mean.nc
spectralcov.nc
)
foreach( file ${saber_test_files} )
    set( fullpath ${CMAKE_CURRENT_SOURCE_DIR}/saber/test/testdata/${file} )
    if( NOT EXISTS ${fullpath} )
        message( STATUS " - Downloading: ${fullpath}" )
        file( DOWNLOAD https://media.githubusercontent.com/media/jcsda/saber/develop/test/testdata/${file} ${fullpath} )
    endif()
endforeach()
set( fullpath ${CMAKE_CURRENT_SOURCE_DIR}/saber/quench/data/landsea.nc )
if( NOT EXISTS ${fullpath} )
    message( STATUS " - Downloading: ${fullpath}" )
    file( DOWNLOAD https://media.githubusercontent.com/media/jcsda/saber/develop/quench/data/landsea.nc ${fullpath} )
endif()

# ECSABER compile definition required for OOPS
if( ECSABER )
  add_compile_definitions( ECSABER )
endif()

# ENABLE_EVIL compile definition required for OOPS
if( ECEVIL )
  add_compile_definitions( ECEVIL )
endif()

# Do not compile quench if the environment variable is not set
if( NOT DEFINED ENABLE_QUENCH )
    set( ENABLE_QUENCH OFF CACHE BOOL "Compile quench" )
endif()

# Add subdirectories
add_subdirectory( oops_patch )
if( ECSABER )
    add_subdirectory( vader )
    add_subdirectory( saber )
endif()
