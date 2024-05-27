/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "quench/Variables.h"

#include "eckit/config/Configuration.h"

#include "oops/util/ConfigFunctions.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
  : vars_(other.vars_) {
  oops::Log::trace() << classname() << "::Variables" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::Variables starting" << std::endl;

  if (util::isVector(config)) {
    vars_ = oops::JediVariables(config.getStringVector("."));
  } else if (config.has("variables")) {
    vars_ = oops::JediVariables(config, "variables");
  } else if (config.has("variables list")) {
    vars_ = oops::JediVariables(config, "variables list");
  } else {
    throw eckit::Exception("wrong variables configuration", Here());
  }

  oops::Log::trace() << classname() << "::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars) {
  oops::Log::trace() << classname() << "::Variables starting" << std::endl;

  vars_ = oops::JediVariables(vars);

  oops::Log::trace() << classname() << "::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << vars_ << std::endl;

  oops::Log::trace() << classname() << "::print end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
