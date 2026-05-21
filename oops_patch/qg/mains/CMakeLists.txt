if( ECEVIL )
    ecbuild_add_executable( TARGET   qg_assim_ensemble_patched.x
                            SOURCES  qgAssimEnsemblePatched.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )

    ecbuild_add_executable( TARGET   qg_bgos.x
                            SOURCES  qgBGOS.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )

    ecbuild_add_executable( TARGET   qg_ensemble_variance.x
                            SOURCES  qgEnsembleVariance.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )

    ecbuild_add_executable( TARGET   qg_4dvar_evil.x
                            SOURCES  qg4DVarEvil.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )

    ecbuild_add_executable( TARGET   qg_make_obs_patched.x
                            SOURCES  qgMakeObsPatched.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )

    ecbuild_add_executable( TARGET   qg_synthetic_data.x
                            SOURCES  qgSyntheticData.cc
                            INCLUDES ${CMAKE_SOURCE_DIR}/oops/qg
                            LIBS    qg
                                    oops_patch
                          )
endif()
