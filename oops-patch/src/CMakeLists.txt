#OOPS
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_jedi_src_files
oops/base/DataSetBase.h
oops/base/FieldSet3D.cc
oops/base/FieldSet3D.h
oops/base/FieldSet4D.cc
oops/base/FieldSet4D.h
oops/base/FieldSets.cc
oops/base/FieldSets.h
oops/base/GeometryData.cc
oops/base/GeometryData.h
oops/base/ModelSpaceCovarianceParametersBase.h
oops/base/ParameterTraitsVariables.cc
oops/base/ParameterTraitsVariables.h
oops/base/Variable.cc
oops/base/Variable.h
oops/base/Variables.cc
oops/base/Variables.h
oops/base/WriteParametersBase.h
oops/generic/AtlasInterpolator.cc
oops/generic/AtlasInterpolator.h
oops/generic/Diffusion.cc
oops/generic/Diffusion.h
oops/generic/gc99.cc
oops/generic/gc99.h
oops/generic/GlobalInterpolator.cc
oops/generic/GlobalInterpolator.h
oops/generic/LocalInterpolatorBase.cc
oops/generic/LocalInterpolatorBase.h
oops/generic/UnstructuredInterpolator.cc
oops/generic/UnstructuredInterpolator.h
oops/interface/ModelData.h
oops/mpi/mpi.cc
oops/mpi/mpi.h
oops/runs/ConvertState.h
oops/util/AnyOf.h
oops/util/AssociativeContainers.h
oops/util/AtlasArrayUtil.cc
oops/util/AtlasArrayUtil.h
oops/util/CompositePath.cc
oops/util/CompositePath.h
oops/util/ConfigFunctions.cc
oops/util/ConfigFunctions.h
oops/util/ConfigHelpers.cc
oops/util/ConfigHelpers.h
oops/util/dateFunctions.cc
oops/util/dateFunctions.h
oops/util/FieldSetHelpers.cc
oops/util/FieldSetHelpers.h
oops/util/FieldSetOperations.cc
oops/util/FieldSetOperations.h
oops/util/FloatCompare.h
oops/util/FunctionSpaceHelpers.cc
oops/util/FunctionSpaceHelpers.h
oops/util/gatherPrint.h
oops/util/IntSetParser.cc
oops/util/IntSetParser.h
oops/util/missingValues.cc
oops/util/missingValues.h
oops/util/NamedEnumerator.h
oops/util/parameters/ConfigurationParameter.cc
oops/util/parameters/ConfigurationParameter.h
oops/util/parameters/EmptyParameters.h
oops/util/parameters/GenericParameters.h
oops/util/parameters/HasDiracParameters_.h
oops/util/parameters/HasParameters_.h
oops/util/parameters/HasReadParameters_.h
oops/util/parameters/HasWriteParameters_.h
oops/util/parameters/IgnoreOtherParameters.cc
oops/util/parameters/IgnoreOtherParameters.h
oops/util/parameters/NumericConstraints.h
oops/util/parameters/ObjectJsonSchema.cc
oops/util/parameters/ObjectJsonSchema.h
oops/util/parameters/OptionalParameter.h
oops/util/parameters/OptionalParameter.cc
oops/util/parameters/OptionalPolymorphicParameter.h
oops/util/parameters/Parameter.h
oops/util/parameters/ParameterBase.cc
oops/util/parameters/ParameterBase.h
oops/util/parameters/ParameterConstraint.h
oops/util/parameters/Parameters.cc
oops/util/parameters/Parameters.h
oops/util/parameters/ParametersOrConfiguration.h
oops/util/parameters/ParameterTraits.cc
oops/util/parameters/ParameterTraits.h
oops/util/parameters/ParameterTraitsAnyOf.h
oops/util/parameters/ParameterTraitsScalarOrMap.h
oops/util/parameters/PolymorphicParameter.h
oops/util/parameters/PolymorphicParameterTraits.h
oops/util/parameters/PropertyJsonSchema.cc
oops/util/parameters/PropertyJsonSchema.h
oops/util/parameters/RequiredParameter.h
oops/util/parameters/RequiredPolymorphicParameter.h
oops/util/PartialDateTime.cc
oops/util/PartialDateTime.h
oops/util/Random.h
oops/util/RandomField.h
oops/util/ScalarOrMap.h
oops/util/Serializable.h
oops/util/stringFunctions.cc
oops/util/stringFunctions.h
oops/util/TypeTraits.h
)

# Copy and update source from OOPS
patch_init( "${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh" )
foreach( src_file ${oops_jedi_src_files} )
    check_diff( ${CMAKE_CURRENT_SOURCE_DIR}/../../oops-jedi/src/${src_file}
                ${CMAKE_CURRENT_SOURCE_DIR}/${src_file}
                ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )

    list( APPEND oops-patch_src_files ${src_file} )
endforeach()
patch_apply( "${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh" )

# Check update summary
add_custom_command( OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/update_result.sh"
                    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/update_patch.sh" )

list( APPEND oops_ecmwf_src_files
util/abor1_cpp.h
util/DateTime.cc
util/DateTime.h
util/Duration.cc
util/Duration.h
util/formats.h
util/kinds.F90
util/Logger.h
util/ObjectCounter.h
util/Printable.h
util/Timer.h
)

foreach( src_file ${oops_ecmwf_src_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/oops/src/${src_file}
                     ${CMAKE_CURRENT_SOURCE_DIR}/oops/${src_file} )

    list( APPEND oops-patch_src_files oops/${src_file} )
endforeach()

list( APPEND util-patch_src_files
util/LibOOPS.h
)

foreach( src_file ${util-patch_src_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/oops/src/${src_file}
                     ${CMAKE_CURRENT_SOURCE_DIR}/${src_file} )
    list( APPEND oops-patch_src_files ${src_file} )
endforeach()

ecbuild_add_library( TARGET          oops-patch
                     SOURCES         ${oops-patch_src_files}
                     PUBLIC_LIBS     NetCDF::NetCDF_Fortran NetCDF::NetCDF_C Eigen3::Eigen eckit fckit atlas_f Boost::boost oops
                     PUBLIC_INCLUDES ${EIGEN3_INCLUDE_DIR}
                                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
                                     $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                                     $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                                     $<INSTALL_INTERFACE:include>
                                     $<INSTALL_INTERFACE:${Fortran_MODULE_INSTALL_DIR}>
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                   )
if( NOT OOPS_LINK_LAPACK )
  # Only here to allow CY47R2.DP to be bitreproducible with operational CY47R1.DP; It can be removed for CY48R1.
  ecbuild_warn( "Not linking oops with LAPACK_LIBRARIES. \nUndefined references need to be resolved later." )
else()
  target_link_libraries( oops PRIVATE ${LAPACK_LIBRARIES} )
endif()
