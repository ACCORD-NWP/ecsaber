/*
 * (C) Copyright 2022 UCAR.
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "mo/eval_dry_air_density.h"
#include "oops/util/Logger.h"
#include "vader/recipes/DryAirDensityLevelsMinusOne.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char DryAirDensityLevelsMinusOne_C::Name[] = "DryAirDensityLevelsMinusOne_C";
const std::vector<std::string> DryAirDensityLevelsMinusOne_C::
    Ingredients = {"height", "height_levels", "potential_temperature",
                   "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
                   "mass_content_of_cloud_ice_in_atmosphere_layer",
                   "specific_humidity", "air_pressure_levels"};
    // For renaming to CCPP standard in (hopefully) near future:
    // Ingredients = {"height_above_mean_sea_level_at_interface",
    //                "height_above_mean_sea_levels_extended_up_by_one",
    //                "air_potential_temperature_at_interface",
    //                "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    //                "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    //                "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    //                "air_pressure_extended_up_by_one"}

// Register the maker
static RecipeMaker<DryAirDensityLevelsMinusOne_C>
    makerDryAirDensityLevelsMinusOne_C_(DryAirDensityLevelsMinusOne_C::Name);

DryAirDensityLevelsMinusOne_C::DryAirDensityLevelsMinusOne_C(const Parameters_ & params,
                                        const VaderConfigVars & configVariables) :
        configVariables_{configVariables}

{
    oops::Log::trace() <<
        "DryAirDensityLevelsMinusOne_C::DryAirDensityLevelsMinusOne_C(params)" << std::endl;
}

std::string DryAirDensityLevelsMinusOne_C::name() const
{
    return DryAirDensityLevelsMinusOne_C::Name;
}

std::string DryAirDensityLevelsMinusOne_C::product() const
{
    return "dry_air_density_levels_minus_one";
    // For renaming to CCPP standard in (hopefully) near future:
    // return "dry_air_density";
}

std::vector<std::string> DryAirDensityLevelsMinusOne_C::ingredients() const
{
    return DryAirDensityLevelsMinusOne_C::Ingredients;
}

size_t DryAirDensityLevelsMinusOne_C::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset["air_pressure_levels"].levels() - 1;
    // For renaming to CCPP standard in (hopefully) near future:
    //  return afieldset["air_pressure_extended_up_by_one"].levels() - 1;
}

atlas::FunctionSpace DryAirDensityLevelsMinusOne_C::productFunctionSpace(const atlas::FieldSet
                                                                 & afieldset) const
{
    return afieldset.field("air_pressure_levels").functionspace();
}

bool DryAirDensityLevelsMinusOne_C::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_C::executeNL function"
      << std::endl;
    mo::eval_new_dry_air_density_from_pressure_levels_nl(afieldset);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_C::executeNL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_C::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_C::executeTL function"
        << std::endl;
    mo::eval_new_dry_air_density_from_pressure_levels_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_C::executeTL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_C::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne::executeAD function"
        << std::endl;
    mo::eval_new_dry_air_density_from_pressure_levels_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_C::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
