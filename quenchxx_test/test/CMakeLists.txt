list( APPEND quenchxx_testinput
  testinput/evil_assim_ensemble.json
  testinput/evil_background.json
  testinput/evil_control_3dvar.json
  testinput/evil_control_update_R.json
  testinput/evil_control_update_S.json
  testinput/evil_dual_3dvar.json
  testinput/evil_dual_update_D.json
  testinput/evil_dual_update_R.json
  testinput/evil_dual_update_S.json
  testinput/evil_ensemble.json
  testinput/evil_ensemble_variance.json
  testinput/evil_makeobs.json
  testinput/evil_primal_3dvar.json
  testinput/evil_primal_update_D.json
  testinput/evil_primal_update_R.json
  testinput/evil_primal_update_S.json
  testinput/evil_stddev.json
  testinput/evil_truth.json
  testinput/glb_3dvar.json
  testinput/glb_3densvar.json
  testinput/glb_4densvar.json
  testinput/glb_background_06.json
  testinput/glb_background_12.json
  testinput/glb_background_18.json
  testinput/glb_ensemble_06.json
  testinput/glb_ensemble_12.json
  testinput/glb_ensemble_18.json
  testinput/glb_hofx.json
  testinput/glb_interfaces.json
  testinput/glb_makeobs_06.json
  testinput/glb_makeobs_12.json
  testinput/glb_makeobs_18.json
  testinput/glb_stddev.json
  testinput/glb_truth_06.json
  testinput/glb_truth_12.json
  testinput/glb_truth_18.json
  testinput/reg_3dvar.json
  testinput/reg_3densvar.json
  testinput/reg_4densvar.json
  testinput/reg_background_06.json
  testinput/reg_background_12.json
  testinput/reg_background_18.json
  testinput/reg_ensemble_06.json
  testinput/reg_ensemble_12.json
  testinput/reg_ensemble_18.json
  testinput/reg_hofx.json
  testinput/reg_interfaces.json
  testinput/reg_makeobs_06.json
  testinput/reg_makeobs_12.json
  testinput/reg_makeobs_18.json
  testinput/reg_stddev.json
  testinput/reg_truth_06.json
  testinput/reg_truth_12.json
  testinput/reg_truth_18.json
)

list( APPEND quenchxx_testref
  testref/evil_assim_ensemble.ref
  testref/evil_background.ref
  testref/evil_control_3dvar.ref
  testref/evil_control_update_R.ref
  testref/evil_control_update_S.ref
  testref/evil_dual_3dvar.ref
  testref/evil_dual_update_D.ref
  testref/evil_dual_update_R.ref
  testref/evil_dual_update_S.ref
  testref/evil_ensemble.ref
  testref/evil_ensemble_variance.ref
  testref/evil_makeobs.ref
  testref/evil_primal_3dvar.ref
  testref/evil_primal_update_D.ref
  testref/evil_primal_update_R.ref
  testref/evil_primal_update_S.ref
  testref/evil_stddev.ref
  testref/evil_truth.ref
  testref/glb_3dvar.ref
  testref/glb_3densvar.ref
  testref/glb_4densvar.ref
  testref/glb_background_06.ref
  testref/glb_background_12.ref
  testref/glb_background_18.ref
  testref/glb_ensemble_06.ref
  testref/glb_ensemble_12.ref
  testref/glb_ensemble_18.ref
  testref/glb_hofx.ref
  testref/glb_makeobs_06.ref
  testref/glb_makeobs_12.ref
  testref/glb_makeobs_18.ref
  testref/glb_stddev.ref
  testref/glb_truth_06.ref
  testref/glb_truth_12.ref
  testref/glb_truth_18.ref
  testref/glb_3dvar.ref
  testref/reg_3densvar.ref
  testref/reg_4densvar.ref
  testref/reg_background_06.ref
  testref/reg_background_12.ref
  testref/reg_background_18.ref
  testref/reg_ensemble_06.ref
  testref/reg_ensemble_12.ref
  testref/reg_ensemble_18.ref
  testref/reg_hofx.ref
  testref/reg_makeobs_06.ref
  testref/reg_makeobs_12.ref
  testref/reg_makeobs_18.ref
  testref/reg_stddev.ref
  testref/reg_truth_06.ref
  testref/reg_truth_12.ref
  testref/reg_truth_18.ref
)

# Create directory for data
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata )

# Create directory for test input and symlink all files
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput )
foreach( FILENAME ${quenchxx_testinput} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
                     ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

# Create directory for reference output and symlink all files
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref )
foreach( FILENAME ${quenchxx_testref} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
                     ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

# Show in project
ecbuild_add_resources( TARGET       quenchxx_test_scripts
                       SOURCES_PACK ${quenchxx_testinput} )

#####################################################################

# Tests coding norms

ecbuild_add_test( TARGET quenchxx_coding_norms_test
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}/interface # TODO: should be extended to interface and internal
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# MPI/OMP setup
set( omp 1 )
foreach( mpi "1" "4" )
#    foreach( DOMAIN "glb" "reg" )
    foreach( DOMAIN "glb" )
        # Tests that create data other tests might use

        ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_stddev_${mpi}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                               ${mpi}
                               ${omp}
                               "${CMAKE_BINARY_DIR}/bin/saber_quench_convertstate.x testinput/${DOMAIN}_stddev.json"
                               "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_stddev.ref"
                          DEPENDS saber_quench_convertstate.x )
        
        foreach( CYCLE "06" "12" "18" )
            ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_truth_${CYCLE}_${mpi}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                              ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                                   ${mpi}
                                   ${omp}
                                   "${CMAKE_BINARY_DIR}/bin/saber_quench_convertstate.x testinput/${DOMAIN}_truth_${CYCLE}.json"
                                   "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_truth_${CYCLE}.ref"
                              DEPENDS saber_quench_convertstate.x )
        
            ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_background_${CYCLE}_${mpi}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                              ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                                   ${mpi}
                                   ${omp}
                                   "${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x testinput/${DOMAIN}_background_${CYCLE}.json"
                                   "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_background_${CYCLE}.ref"
                              DEPENDS saber_quench_error_covariance_toolbox.x )
        
            ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_ensemble_${CYCLE}_${mpi}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                              ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                                   ${mpi}
                                   ${omp}
                                   "${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x testinput/${DOMAIN}_ensemble_${CYCLE}.json"
                                   "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_ensemble_${CYCLE}.ref"
                              DEPENDS saber_quench_error_covariance_toolbox.x )
        
            ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_makeobs_${CYCLE}_${mpi}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                              ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                                   ${mpi}
                                   ${omp}
                                   "${CMAKE_BINARY_DIR}/bin/saber_quench_makeobs.x testinput/${DOMAIN}_makeobs_${CYCLE}.json"
                                   "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_makeobs_${CYCLE}.ref"
                              DEPENDS saber_quench_makeobs.x )
        endforeach()
    
        # Test internal classes
    
#        ecbuild_add_test( TARGET quench_internal_resolution
#                          SOURCES quench/Resolution.cc quench/TestConfig.h
#                          ARGS "testinput/${DOMAIN}_test.json"
#                          LIBS    quench eckit)
    
        # Test interface classes with quench
        
        if( OFF )
            ecbuild_add_test( TARGET  quenchxx_interface_error_covariance_${mpi}
                              SOURCES interface/TestErrorCovariance.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_geometry_${mpi}
                              SOURCES interface/TestGeometry.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_increment_${mpi}
                              SOURCES interface/TestIncrement.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_linear_model_${mpi}
                              SOURCES interface/TestLinearModel.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_localization_${mpi}
                              SOURCES interface/TestLocalization.cc
                              ARGS   "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_locations_${mpi}
                              SOURCES interface/TestLocations.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_model_${mpi}
                              SOURCES interface/TestModel.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_model_aux_control_${mpi}
                              SOURCES interface/TestModelAuxControl.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_model_aux_covariance_${mpi}
                              SOURCES interface/TestModelAuxCovariance.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
    #        ecbuild_add_test( TARGET  quenchxx_interface_model_aux_increment
    #                          SOURCES interface/TestModelAuxIncrement.cc
    #                          ARGS    "testinput/${DOMAIN}_interfaces.json"
    #                          LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_obs_error_covariance_${mpi}
                              SOURCES interface/TestObsErrorCovariance.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_observation_space_${mpi}
                              SOURCES interface/TestObservationSpace.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_obs_vector_${mpi}
                              SOURCES interface/TestObsVector.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_state_${mpi}
                              SOURCES interface/TestState.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
            
            ecbuild_add_test( TARGET  quenchxx_interface_variables_${mpi}
                              SOURCES interface/TestVariables.cc
                              ARGS    "testinput/${DOMAIN}_interfaces.json"
                              LIBS    quench )
        endif()
        
        # Test oops applications with quench
        
        ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_hofx_${mpi}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                               ${mpi}
                               ${omp}
                               "${CMAKE_BINARY_DIR}/bin/saber_quench_hofx.x testinput/${DOMAIN}_hofx.json"
                               "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_hofx.ref"
                          DEPENDS saber_quench_hofx.x )
        
        ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_3dvar_${mpi}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                               ${mpi}
                               ${omp}
                               "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/${DOMAIN}_3dvar.json"
                               "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_3dvar.ref"
                          DEPENDS saber_quench_variational.x )
        
        ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_3densvar_${mpi}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                               ${mpi}
                               ${omp}
                               "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/${DOMAIN}_3densvar.json"
                               "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_3densvar.ref"
                          DEPENDS saber_quench_variational.x )
        
        ecbuild_add_test( TARGET quenchxx_test_${DOMAIN}_4densvar_${mpi}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                          ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                               ${mpi}
                               ${omp}
                               "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/${DOMAIN}_4densvar.json"
                               "${CMAKE_CURRENT_BINARY_DIR}/testref/${DOMAIN}_4densvar.ref"
                          DEPENDS saber_quench_variational.x )
    endforeach()

    # EVIL tests  
    ecbuild_add_test( TARGET quenchxx_test_evil_truth_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_convertstate.x testinput/evil_truth.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_truth.ref"
                      DEPENDS saber_quench_convertstate.x )
    
    ecbuild_add_test( TARGET quenchxx_test_evil_stddev_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_convertstate.x testinput/evil_stddev.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_stddev.ref"
                      DEPENDS saber_quench_convertstate.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_background_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x testinput/evil_background.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_background.ref"
                      DEPENDS saber_quench_error_covariance_toolbox.x )
    
    ecbuild_add_test( TARGET quenchxx_test_evil_ensemble_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x testinput/evil_ensemble.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_ensemble.ref"
                      DEPENDS saber_quench_error_covariance_toolbox.x )
    
    ecbuild_add_test( TARGET quenchxx_test_evil_makeobs_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_makeobs.x testinput/evil_makeobs.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_makeobs.ref"
                      DEPENDS saber_quench_makeobs.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_control_3dvar_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/evil_control_3dvar.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_control_3dvar.ref"
                      DEPENDS saber_quench_variational.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_control_update_R_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_control_update_R.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_control_update_R.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_control_update_S_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_control_update_S.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_control_update_S.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_dual_3dvar_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/evil_dual_3dvar.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_dual_3dvar.ref"
                      DEPENDS saber_quench_variational.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_dual_update_D_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_dual_update_D.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_dual_update_D.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_dual_update_R_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_dual_update_R.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_dual_update_R.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_dual_update_S_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_dual_update_S.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_dual_update_S.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_primal_3dvar_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_variational.x testinput/evil_primal_3dvar.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_primal_3dvar.ref"
                      DEPENDS saber_quench_variational.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_primal_update_D_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_primal_update_D.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_primal_update_D.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_primal_update_R_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_primal_update_R.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_primal_update_R.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_primal_update_S_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_evil_update.x testinput/evil_primal_update_S.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_primal_update_S.ref"
                      DEPENDS saber_quench_evil_update.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_assim_ensemble_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_assim_ensemble.x testinput/evil_assim_ensemble.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_assim_ensemble.ref"
                      DEPENDS saber_quench_assim_ensemble.x )

    ecbuild_add_test( TARGET quenchxx_test_evil_ensemble_variance_${mpi}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.sh
                      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/ecmwf_compare.py
                           ${mpi}
                           ${omp}
                           "${CMAKE_BINARY_DIR}/bin/saber_quench_ensemble_variance.x testinput/evil_ensemble_variance.json"
                           "${CMAKE_CURRENT_BINARY_DIR}/testref/evil_ensemble_variance.ref"
                      DEPENDS saber_quench_ensemble_variance.x )
endforeach()
