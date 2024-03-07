/*
 * (C) Copyright 2024- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>
#include <vector>

#include "atlas/field/Field.h"
#include "mo/eval_air_pressure_levels.h"
#include "oops/util/Logger.h"
#include "vader/recipes/AirPressureExtendedUpByOne.h"

namespace vader {

// -------------------------------------------------------------------------------------------------
// Static attribute initialization
const char AirPressureExtendedUpByOne_A::Name[] = "AirPressureExtendedUpByOne_A";

// Register the maker
static RecipeMaker<AirPressureExtendedUpByOne_A>
  makerAirPressureExtendedUpByOne_A_(AirPressureExtendedUpByOne_A::Name);

// -------------------------------------------------------------------------------------------------

AirPressureExtendedUpByOne_A::AirPressureExtendedUpByOne_A(const Parameters_ & params,
                                         const VaderConfigVars &)
{
    oops::Log::trace() << "AirPressureExtendedUpByOne_A::AirPressureExtendedUpByOne_A Done"
                       << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> AirPressureExtendedUpByOne_A::ingredients() const {
    return {"exner_levels_minus_one", "air_pressure_levels_minus_one",
            "potential_temperature", "height_levels"};
    // CCPP Names for (hopefully near) future
    // return {"dimensionless_exner_function", "air_pressure",
    //         "air_potential_temperature_at_interface",
    //         "height_above_mean_sea_level_extended_up_by_one"};
}

// -------------------------------------------------------------------------------------------------

size_t AirPressureExtendedUpByOne_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels_minus_one").shape(1) + 1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressureExtendedUpByOne_A::productFunctionSpace(
  const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels_minus_one").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressureExtendedUpByOne_A::executeNL(atlas::FieldSet & afieldset) {
    oops::Log::trace() << "AirPressureExtendedUpByOne_A::executeNL Starting" << std::endl;
    mo::eval_air_pressure_levels_nl(afieldset);
    oops::Log::trace() << "AirPressureExtendedUpByOne_A::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
