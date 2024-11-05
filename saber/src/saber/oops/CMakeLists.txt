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

# Utilities
Utilities.h
Utilities.cc

# Factories
instantiateCovarFactory.h
instantiateLocalizationFactory.h

# Utilities
ECUtilities.cc
ECUtilities.h

# Other applications
AssimEnsemble.h
EnsembleVariance.h
MakeObs.h
BGOS.h
SyntheticData.h

# EVIL
DRPLanczosEVILMinimizer.h
EvilUpdate.h
instantiateMinFactory.h
RitzPairs.h
PLanczosEVILMinimizer.h
RPLanczosEVILMinimizer.h
SQRTPLanczosEVILMinimizer.h

# Temporary versions for cy46
DRPLanczosEVILMinimizer_cy46.h
EvilUpdate_cy46.h
instantiateMinFactory_cy46.h
RitzPairs_cy46.h
)

set( oops_src_files ${oops_src_files_list} PARENT_SCOPE )
