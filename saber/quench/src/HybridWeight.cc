/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/HybridWeight.h"

#include <cmath>
#include <vector>

#include "util/Logger.h"

#include "src/Increment.h"

namespace quench {

// -----------------------------------------------------------------------------

HybridWeight::HybridWeight(const eckit::Configuration & config)
  : wgt_(std::sqrt(config.getDouble("weight"))) {
  oops::Log::trace() << "HybridWeight::HybridWeight done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
