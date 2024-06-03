# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

# Update flag
set( ECSABER_UPDATE ON )
if( DEFINED ENV{ECSABER_UPDATE} )
    set( ECSABER_UPDATE $ENV{ECSABER_UPDATE} )
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

function( check_diff projectName srcPath dstPath commandPath )
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
set( ENABLE_QUENCH ON )
if( DEFINED ENV{ENABLE_QUENCH} )
    set( ENABLE_QUENCH $ENV{ENABLE_QUENCH} )
endif()

# Add subdirectories
if( ECSABER )
    add_subdirectory( oops-patch )
    add_subdirectory( vader )
    add_subdirectory( saber )
    if( ENABLE_QUENCH )
        add_subdirectory( quench )
    endif()
endif()
