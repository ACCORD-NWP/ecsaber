/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/HybridWeight.h"

#include <cmath>
#include <vector>

#include "src/Increment.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

HybridWeight::HybridWeight(const eckit::Configuration & config)
  : wgt_(std::sqrt(config.getDouble("weight"))) {
  oops::Log::trace() << classname() << "::HybridWeight" << std::endl;
}


// -----------------------------------------------------------------------------

void HybridWeight::multiply(Increment & dx) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  dx *= wgt_;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
