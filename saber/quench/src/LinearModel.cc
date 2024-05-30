/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Increment.h"
#include "src/LinearModel.h"
#include "src/TraitsFwd.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

static oops::LinearModelMaker<Traits, LinearModel> makerLinearModelDefault_("default");

// -----------------------------------------------------------------------------

LinearModel::LinearModel(const Geometry &,
                         const eckit::Configuration & config)
  : timeResolution_(config.getString("tstep", "PT6H")) {
  oops::Log::trace() << classname() << "::LinearModel starting" << std::endl;

  oops::Log::info() << "Persistance linear model" << std::endl;

  oops::Log::trace() << classname() << "::LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearModel::stepTL(Increment & dx,
                         const ModelAuxIncrement &) const {
  oops::Log::trace() << classname() << "::stepTL starting" << std::endl;

  dx.validTime() += timeResolution_;

  oops::Log::trace() << classname() << "::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearModel::stepAD(Increment & dx,
                         ModelAuxIncrement &) const {
  oops::Log::trace() << classname() << "::stepAD starting" << std::endl;

  dx.validTime() -= timeResolution_;

  oops::Log::trace() << classname() << "::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
