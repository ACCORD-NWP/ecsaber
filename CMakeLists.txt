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
function( check_diff projectName srcPath dstPath commandPath )
    if( ECSABER_UPDATE )
        execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff.sh ${srcPath} ${dstPath} ${commandPath} )
    endif()
endfunction()

function( check_diff_recursively input commandPath )
    if( ECSABER_UPDATE )
        execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff_recursively.sh ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${input} ${commandPath} )
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
    if( ENABLE_QUENCH )
        add_subdirectory( quench )
    endif()
    add_subdirectory( saber )
endif()
