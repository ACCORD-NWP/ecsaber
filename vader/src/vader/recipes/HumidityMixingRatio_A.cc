/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/for_each.h"
#include "atlas/util/Metadata.h"
#include "oops/util/Logger.h"
#include "vader/recipes/HumidityMixingRatio.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char HumidityMixingRatio_A::Name[] = "HumidityMixingRatio_A";
const oops::JediVariables HumidityMixingRatio_A::Ingredients{
                      {"water_vapor_mixing_ratio_wrt_moist_air"}};

// Register the maker
static RecipeMaker<HumidityMixingRatio_A> makerHumidityMixingRatio_(HumidityMixingRatio_A::Name);

HumidityMixingRatio_A::HumidityMixingRatio_A(const Parameters_ & params,
                                       const VaderConfigVars & configVariables) :
    configVariables_{configVariables}
{
    oops::Log::trace() << "HumidityMixingRatio_A::HumidityMixingRatio_A(params)" << std::endl;
}

std::string HumidityMixingRatio_A::name() const
{
    return HumidityMixingRatio_A::Name;
}

oops::Variable HumidityMixingRatio_A::product() const
{
    return oops::Variable{"humidity_mixing_ratio"};
}

oops::JediVariables HumidityMixingRatio_A::ingredients() const
{
    return HumidityMixingRatio_A::Ingredients;
}

size_t HumidityMixingRatio_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("water_vapor_mixing_ratio_wrt_moist_air").shape(1);
}

atlas::FunctionSpace HumidityMixingRatio_A::productFunctionSpace
                                              (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("water_vapor_mixing_ratio_wrt_moist_air").functionspace();
}

bool HumidityMixingRatio_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering HumidityMixingRatio_A::executeNL function" << std::endl;

    atlas::field::for_each_value(afieldset["water_vapor_mixing_ratio_wrt_moist_air"],
                                 afieldset["humidity_mixing_ratio"],
                                 [&](const double q, double& mixr) {
        mixr = q / (1. - q);
    });

    oops::Log::trace() << "leaving HumidityMixingRatio_A::executeNL function" << std::endl;

    return true;
}

}  // namespace vader
