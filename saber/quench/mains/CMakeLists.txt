# (C) Copyright 2017-2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Executables
ecbuild_add_executable( TARGET  saber_quench_error_covariance_toolbox.x
                        SOURCES quenchErrorCovarianceToolbox.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_process_perts.x
                        SOURCES quenchProcessPerts.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                        SOURCES quenchConvertState.cc
                        LIBS    quench
                                vader
                                saber )

ecbuild_add_executable( TARGET  saber_quench_hofx.x
                        SOURCES quenchHofX.cc
                        LIBS    quench )

ecbuild_add_executable( TARGET  saber_quench_makeobs.x
                        SOURCES quenchMakeObs.cc
                        LIBS    quench )

ecbuild_add_executable( TARGET  saber_quench_variational.x
                        SOURCES quenchVariational.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_evil_update.x
                        SOURCES quenchEvilUpdate.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_assim_ensemble.x
                        SOURCES quenchAssimEnsemble.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_ensemble_variance.x
                        SOURCES quenchEnsembleVariance.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_bgos.x
                        SOURCES quenchBGOS.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )
