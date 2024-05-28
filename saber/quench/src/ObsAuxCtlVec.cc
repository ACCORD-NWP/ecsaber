/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ObsAuxCtlVec.h"

#include "oops/interface/ObsAuxCovarianceBase.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

static oops::ObsAuxCtlVecMaker<Traits, ObsAuxCtlVec> makeObsAuxCtlVecDefault_("default");

// -----------------------------------------------------------------------------

}  // namespace quench
