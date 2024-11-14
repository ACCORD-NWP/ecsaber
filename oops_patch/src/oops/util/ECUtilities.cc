/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ECUtilities.h"

#include "oops/util/ConfigFunctions.h"

namespace util {

// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration & conf,
               const size_t & member) {
  oops::Log::trace() << "setMember starting" << std::endl;

  if (conf.has("member pattern")) {
    std::string memberPattern = conf.getString("member pattern");
    seekAndReplace(conf, memberPattern, std::to_string(member));
  } else {
    conf.set("member", member);
  }

  oops::Log::trace() << "setMember done" << std::endl;
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration templatedVarsConf(const oops::JediVariables & vars) {
  eckit::LocalConfiguration varConf;
  varConf.set("variables list", vars.variables());
  return varConf;
}

// -----------------------------------------------------------------------------

}  // namespace util
