# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

# Update flag
set( ECSABER_UPDATE OFF )
if( DEFINED ENV{ECSABER_UPDATE} )
    set( ECSABER_UPDATE $ENV{ECSABER_UPDATE} )
endif()

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
function( patch_init commandPath )
    if( ECSABER_UPDATE )
        # Remove update script
        file( REMOVE ${commandPath} )

        # Create update script
        file( WRITE ${commandPath} "#!/usr/bin/env bash\n" )
    endif()
endfunction()

function( check_diff srcPath dstPath commandPath )
    if( ECSABER_UPDATE )
        # Set JEDIPREFIX flag
        set( JEDIPREFIX ON )

        # Run check_diff script
        execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff.sh ${srcPath} ${dstPath} ${commandPath} ${JEDIPREFIX} )
    endif()
endfunction()

function( check_diff_recursively input commandPath )
    if( ECSABER_UPDATE )
        # Set JEDIPREFIX flag
        set( EXTRA_ARGS ${ARGN})    
        list( LENGTH EXTRA_ARGS EXTRA_COUNT )
        if( ${EXTRA_COUNT} GREATER 0)
            list( GET EXTRA_ARGS 0 JEDIPREFIX )
        else()
            set( JEDIPREFIX ON )
        endif()

        # Run check_diff_recursively script
        execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff_recursively.sh ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${input} ${commandPath} ${JEDIPREFIX} )
    endif()
endfunction()

function( patch_apply commandPath )
    if( ECSABER_UPDATE )
        # Run update script
        execute_process( COMMAND bash ${commandPath}
                         COMMAND_ERROR_IS_FATAL ANY )

        # Remove update script
        file( REMOVE ${commandPath} )
    endif()
endfunction()

# QUENCH compilation
set( ENABLE_QUENCHXX ON )
if( DEFINED ENV{ENABLE_QUENCHXX} )
    set( ENABLE_QUENCHXX $ENV{ENABLE_QUENCHXX} )
endif()

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
ptheta_bar_mean.nc
Prho_bar_Mean.nc
spectralcov.nc
)
foreach( file ${saber_test_files} )
    set( fullpath ${CMAKE_SOURCE_DIR}/ecsaber/saber/test/testdata/${file} )
    if( NOT EXISTS ${fullpath} )
        message( STATUS " - Downloading: ${fullpath}" )
        file( DOWNLOAD https://media.githubusercontent.com/media/jcsda/saber/develop/test/testdata/${file} ${fullpath} )
    endif()
endforeach()
set( fullpath ${CMAKE_SOURCE_DIR}/ecsaber/saber/quench/data/landsea.nc )
if( NOT EXISTS ${fullpath} )
    message( STATUS " - Downloading: ${fullpath}" )
    file( DOWNLOAD https://media.githubusercontent.com/media/jcsda/saber/develop/quench/data/landsea.nc ${fullpath} )
endif()

# Add subdirectories
add_subdirectory( oops-patch )
add_subdirectory( vader )
add_subdirectory( saber )
if( ENABLE_QUENCHXX )
    add_subdirectory( quenchxx_test )
endif()
