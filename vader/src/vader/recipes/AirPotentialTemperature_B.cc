/*
 * (C) Copyright 2023  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "oops/util/Logger.h"
#include "vader/recipes/AirPotentialTemperature.h"

namespace vader
{
// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPotentialTemperature_B::Name[] = "AirPotentialTemperature_B";
const oops::JediVariables AirPotentialTemperature_B::Ingredients{
      std::vector<std::string>{"air_temperature", "air_pressure_to_kappa"}};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPotentialTemperature_B> makerTempToPTemp_(AirPotentialTemperature_B::Name);

// -------------------------------------------------------------------------------------------------

AirPotentialTemperature_B::AirPotentialTemperature_B(const Parameters_ & params,
                                        const VaderConfigVars & configVariables) {
    oops::Log::trace() << "AirPotentialTemperature_B::AirPotentialTemperature_B(params)"
        << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPotentialTemperature_B::name() const {
    return AirPotentialTemperature_B::Name;
}

// -------------------------------------------------------------------------------------------------

oops::Variable AirPotentialTemperature_B::product() const {
    return oops::Variable("potential_temperature");
}

// -------------------------------------------------------------------------------------------------

oops::JediVariables AirPotentialTemperature_B::ingredients() const {
    return AirPotentialTemperature_B::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPotentialTemperature_B::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_temperature").shape(1);
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPotentialTemperature_B::productFunctionSpace
                                                (const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_temperature").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPotentialTemperature_B::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering AirPotentialTemperature_B::executeNL function" << std::endl;

    // Get fields
    atlas::Field air_temperature_field = afieldset.field("air_temperature");
    atlas::Field air_pressure_to_kappa_field = afieldset.field("air_pressure_to_kappa");
    atlas::Field potential_temperature_field = afieldset.field("potential_temperature");

    // Get field arrays
    auto air_temperature = atlas::array::make_view<double, 2>(air_temperature_field);
    auto air_pressure_to_kappa = atlas::array::make_view<double, 2>(air_pressure_to_kappa_field);
    auto potential_temperature = atlas::array::make_view<double, 2>(potential_temperature_field);

    // Grid dimensions
    size_t h_size = air_temperature_field.shape(0);
    int v_size = air_temperature_field.shape(1);

    // Calculate the output variable
    for (int vv = 0; vv < v_size; ++vv) {
      for ( size_t hh = 0; hh < h_size ; ++hh ) {
        potential_temperature(hh, vv) = air_temperature(hh, vv) / air_pressure_to_kappa(hh, vv);
      }
    }

    oops::Log::trace() << "leaving AirPotentialTemperature_B::executeNL function" << std::endl;

    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader
