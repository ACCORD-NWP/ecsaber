/*
 * (C) Copyright 2022 UCAR.
 * (C) Crown Copyright 2023-2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "mo/eval_air_density.h"
#include "oops/util/Logger.h"
#include "vader/recipes/AirDensityLevelsMinusOne.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirDensityLevelsMinusOne_B::Name[] = "AirDensityLevelsMinusOne_B";
const oops::JediVariables AirDensityLevelsMinusOne_B::
    Ingredients{std::vector<std::string>{"height", "height_levels", "virtual_potential_temperature",
                                             "air_pressure_levels"}};

// Register the maker
static RecipeMaker<AirDensityLevelsMinusOne_B>
    makerAirDensityLevelsMinusOne_B_(AirDensityLevelsMinusOne_B::Name);

AirDensityLevelsMinusOne_B::AirDensityLevelsMinusOne_B(const Parameters_ & params,
                                        const VaderConfigVars & configVariables) :
        configVariables_{configVariables}

{
    oops::Log::trace() <<
        "AirDensityLevelsMinusOne_B::AirDensityLevelsMinusOne_B(params)" << std::endl;
}

std::string AirDensityLevelsMinusOne_B::name() const
{
    return AirDensityLevelsMinusOne_B::Name;
}

oops::Variable AirDensityLevelsMinusOne_B::product() const
{
    return oops::Variable("air_density_levels_minus_one");
}

oops::JediVariables AirDensityLevelsMinusOne_B::ingredients() const
{
    return AirDensityLevelsMinusOne_B::Ingredients;
}

size_t AirDensityLevelsMinusOne_B::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset["air_pressure_levels"].shape(1) - 1;
}

atlas::FunctionSpace AirDensityLevelsMinusOne_B::productFunctionSpace(const atlas::FieldSet
                                                                 & afieldset) const
{
    return afieldset.field("air_pressure_levels").functionspace();
}

bool AirDensityLevelsMinusOne_B::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering AirDensityLevelsMinusOne_B::executeNL function"
      << std::endl;
    mo::eval_air_density_from_pressure_levels_nl(afieldset);
    oops::Log::trace() << "leaving AirDensityLevelsMinusOne_B::executeNL function" << std::endl;

    return true;
}

bool AirDensityLevelsMinusOne_B::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering AirDensityLevelsMinusOne_B::executeTL function"
        << std::endl;
    mo::eval_air_density_from_pressure_levels_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving AirDensityLevelsMinusOne_B::executeTL function" << std::endl;

    return true;
}

bool AirDensityLevelsMinusOne_B::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering AirDensityLevelsMinusOne::executeAD function"
        << std::endl;
    mo::eval_air_density_from_pressure_levels_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving AirDensityLevelsMinusOne_B::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
