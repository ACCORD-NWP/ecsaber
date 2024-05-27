/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/default/FieldsIODefault.h"

#include <vector>

#include "oops/util/FieldSetHelpers.h"

#include "quench/Geometry.h"
#include "quench/Variables.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

static FieldsIOMaker<FieldsIODefault> makerDefault_("default");

// -----------------------------------------------------------------------------

void FieldsIODefault::read(const Geometry & geom,
                           const Variables & vars,
                           const eckit::Configuration & conf,
                           atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Create variableSizes
  std::vector<size_t> variableSizes;
  for (const auto & var : vars.variablesList()) {
    variableSizes.push_back(geom.levels(var));
  }

  // Read fieldset
  util::readFieldSet(geom.getComm(),
                     geom.functionSpace(),
                     variableSizes,
                     vars.variablesList(),
                     conf,
                     fset);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsIODefault::write(const Geometry & geom,
                            const eckit::Configuration & conf,
                            const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Write fieldset
  util::writeFieldSet(geom.getComm(), conf, fset);

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
