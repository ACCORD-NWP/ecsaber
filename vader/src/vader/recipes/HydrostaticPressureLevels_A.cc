/*
 * (C) Crown Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <iostream>
#include <vector>

#include "atlas/field/Field.h"
#include "mo/control2analysis_varchange.h"
#include "oops/util/Logger.h"
#include "vader/recipes/HydrostaticPressureLevels.h"


namespace vader {


// -------------------------------------------------------------------------------------------------
// Static attribute initialization
const char HydrostaticPressureLevels_A::Name[] =
  "HydrostaticPressureLevels_A";

// Register the maker
static RecipeMaker<HydrostaticPressureLevels_A>
  makerHydrostaticPressureLevels_A_(HydrostaticPressureLevels_A::Name);

// -------------------------------------------------------------------------------------------------

HydrostaticPressureLevels_A::HydrostaticPressureLevels_A(
  const Parameters_ & params,
  const VaderConfigVars & configVariables) {
  oops::Log::trace() <<
    "HydrostaticPressureLevels_A::HydrostaticPressureLevels_A"
    << std::endl;
}

std::string HydrostaticPressureLevels_A::name() const {
  return HydrostaticPressureLevels_A::Name;
}

std::string HydrostaticPressureLevels_A::product() const {
  return "hydrostatic_pressure_levels";
}

std::vector<std::string> HydrostaticPressureLevels_A::ingredients() const {
  return {"hydrostatic_exner_levels"};
}

size_t HydrostaticPressureLevels_A::productLevels(const atlas::FieldSet & afieldset) const {
  return afieldset.field("hydrostatic_exner_levels").shape(1);
}

atlas::FunctionSpace HydrostaticPressureLevels_A::productFunctionSpace(
  const atlas::FieldSet & afieldset) const {
  return afieldset.field("hydrostatic_exner_levels").functionspace();
}

bool HydrostaticPressureLevels_A::executeNL(atlas::FieldSet & afieldset) {
  oops::Log::trace() << "HydrostaticPressureLevels_A::executeNL starting ..." << std::endl;
  mo::evalHydrostaticPressureLevels(afieldset);

  oops::Log::trace() << "HydrostaticPressureLevels_A::executeNL ... done" << std::endl;
  return true;
}


// -------------------------------------------------------------------------------------------------

}  // namespace vader
