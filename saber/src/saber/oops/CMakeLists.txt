# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND oops_src_files_list
# Parameters
ErrorCovarianceParameters.h

# Drivers
ErrorCovariance.h
Localization.h

# Applications
ErrorCovarianceToolbox.h
ProcessPerts.h
AssimEnsemble.h
EnsembleVariance.h
EvilUpdate.h
MakeObs.h
BGOS.h

# Utilities
ECUtilities.cc
ECUtilities.h
Utilities.h
Utilities.cc

# Factories
instantiateCovarFactory.h
instantiateLocalizationFactory.h
instantiateMinFactory.h

# EVIL
RitzPairs.h
PLanczosEVILMinimizer.h
RPLanczosEVILMinimizer.h
SQRTPLanczosEVILMinimizer.h

# Temporary versions for cy46
BGOS_cy46.h
DRPLanczosEVILMinimizer_cy46.h
EvilUpdate_cy46.h
instantiateMinFactory_cy46.h
RitzPairs_cy46.h
SyntheticData_cy46.h
)

set( oops_src_files ${oops_src_files_list} PARENT_SCOPE )
