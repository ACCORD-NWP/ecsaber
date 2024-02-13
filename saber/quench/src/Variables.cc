/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Variables.h"

#include "eckit/config/Configuration.h"

#include "oops/util/ConfigFunctions.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & config) {
  oops::Log::trace() << "Variables::Variables starting" << std::endl;

  if (util::isVector(config)) {
    vars_ = oops::patch::Variables(config.getStringVector("."));
  } else if (config.has("variables")) {
    vars_ = oops::patch::Variables(config, "variables");
  } else if (config.has("variables list")) {
    vars_ = oops::patch::Variables(config, "variables list");
  } else {
    throw eckit::Exception("wrong variables configuration", Here());
  }

  oops::Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
