# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# VADER libraries
if( vader_FOUND )
    set( vader_LIBRARIES "vader" )
else()
    set( vader_LIBRARIES "" )
endif()

# Executables
ecbuild_add_executable( TARGET  saber_quench_error_covariance_toolbox.x
                        SOURCES quenchErrorCovarianceToolbox.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_process_perts.x
                        SOURCES quenchProcessPerts.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                        SOURCES quenchConvertState.cc
                        LIBS    quench )

ecbuild_add_executable( TARGET  saber_quench_hofx.x
                        SOURCES quenchHofX.cc
                        LIBS    quench )

ecbuild_add_executable( TARGET  saber_quench_makeobs.x
                        SOURCES quenchMakeObs.cc
                        LIBS    quench )

ecbuild_add_executable( TARGET  saber_quench_variational.x
                        SOURCES quenchVariational.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )
