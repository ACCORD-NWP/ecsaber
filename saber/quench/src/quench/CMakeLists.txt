# (C) Copyright 2023 Meteorologisk Institutt

list( APPEND quench_src_files

Covariance.h
Fields.cc
Fields.h
FieldsIOBase.cc
FieldsIOBase.h
Geometry.cc
Geometry.h
GeoVaLs.cc
GeoVaLs.h
HorizScaleDecomposition.h
HybridWeight.cc
HybridWeight.h
Increment.cc
Increment.h
Interpolation.cc
Interpolation.h
IncrEnsCtlVec.h
IncrModCtlVec.h
instantiateQuenchMatrices.h
Interpolator.h
LinearModel.cc
LinearModel.h
LinearObsOperator.cc
LinearObsOperator.h
LocalizationMatrix.h
Locations.cc
Locations.h
Logbook.cc
Logbook.h
Model.cc
Model.h
ModelAuxControl.h
ModelAuxCovariance.h
ModelAuxCtlVec.h
ModelAuxControlEstimator.h
ModelAuxIncrement.h
ObsAuxControl.cc
ObsAuxControl.h
ObsAuxCovariance.cc
ObsAuxCovariance.h
ObsAuxCtlVec.cc
ObsAuxCtlVec.h
ObsAuxIncrement.cc
ObsAuxIncrement.h
ObsError.cc
ObsError.h
ObsOperator.cc
ObsOperator.h
ObsSpace.cc
ObsSpace.h
ObsVec.cc
ObsVec.h
State.cc
State.h
Traits.h
Variables.h
Variables.cc

default/FieldsIODefault.cc
default/FieldsIODefault.h

)

ecbuild_add_library( TARGET   quench
                     SOURCES  ${quench_src_files}
                     PUBLIC_LIBS  oops oops-patch saber
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quench
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

target_link_libraries( quench PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )

#Configure include directory layout for build-tree to match install-tree
set(QUENCH_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(quench_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber"
                                     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber/quench)
add_dependencies(quench quench_headers)
target_include_directories(quench PUBLIC $<BUILD_INTERFACE:${QUENCH_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../../..>)
