/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/constants.h"
#include "oops/util/Logger.h"
#include "vader/recipes/AirTemperature.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirTemperature_C::Name[] = "AirTemperature_C";
const char AT[] = "air_temperature";
const char PAPT[] = "perturbation_air_potential_temperature";
const char AP[] = "air_pressure";
const oops::JediVariables AirTemperature_C::Ingredients{std::vector<std::string>{PAPT, AP}};

// Register the maker
static RecipeMaker<AirTemperature_C> makerAirTemperature_C_(AirTemperature_C::Name);

AirTemperature_C::AirTemperature_C(const Parameters_ & params,
                                   const VaderConfigVars & configVariables):
                                   configVariables_{configVariables}
{
    oops::Log::trace() << "AirTemperature_C::AirTemperature_C(params)" << std::endl;
}

std::string AirTemperature_C::name() const
{
    return AirTemperature_C::Name;
}

oops::Variable AirTemperature_C::product() const
{
    return oops::Variable("air_temperature");
}

oops::JediVariables AirTemperature_C::ingredients() const
{
    return AirTemperature_C::Ingredients;
}

size_t AirTemperature_C::productLevels(const atlas::FieldSet & fields) const
{
    return fields[PAPT].levels();
}

atlas::FunctionSpace AirTemperature_C::productFunctionSpace(const atlas::FieldSet & fields) const
{
    return fields[PAPT].functionspace();
}

bool AirTemperature_C::executeNL(atlas::FieldSet & fields)
{
    oops::Log::trace() << "entering AirTemperature_C::executeNL function" << std::endl;

    const double pt_base = configVariables_.getDouble("base_air_potential_temperature");

    auto air_pressure_view = make_view<const double, 2>(fields[AP]);
    auto perturbation_potential_temperature_view = make_view<const double, 2>(fields[PAPT]);
    auto temp_view = make_view<double, 2>(fields[AT]);

    for (idx_t jn = 0; jn < fields[AT].shape(0) ; ++jn) {
      for (idx_t jl = 0; jl < fields[AT].shape(1); ++jl) {
        temp_view(jn, jl) = (perturbation_potential_temperature_view(jn, jl) + pt_base) *
                            pow(air_pressure_view(jn, jl) / mo::constants::p_zero,
                            mo::constants::rd_over_cp);
      }
    }
    oops::Log::trace() << "leaving AirTemperature_C::executeNL function" << std::endl;
    return true;
}

}  // namespace vader
