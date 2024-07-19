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
const char HumidityMixingRatio_B::Name[] = "HumidityMixingRatio_B";
const oops::JediVariables HumidityMixingRatio_B::Ingredients{{"humidity_mixing_ratio_kgkg"}};

// Register the maker
static RecipeMaker<HumidityMixingRatio_B> makerHumidityMixingRatio_(HumidityMixingRatio_B::Name);

HumidityMixingRatio_B::HumidityMixingRatio_B(const Parameters_ & params,
                                       const VaderConfigVars & configVariables) :
    configVariables_{configVariables}
{
    oops::Log::trace() << "HumidityMixingRatio_B::HumidityMixingRatio_B(params)" << std::endl;
}

std::string HumidityMixingRatio_B::name() const
{
    return HumidityMixingRatio_B::Name;
}

oops::Variable HumidityMixingRatio_B::product() const
{
    return oops::Variable{"humidity_mixing_ratio"};
}

oops::JediVariables HumidityMixingRatio_B::ingredients() const
{
    return HumidityMixingRatio_B::Ingredients;
}

size_t HumidityMixingRatio_B::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("humidity_mixing_ratio_kgkg").shape(1);
}

atlas::FunctionSpace HumidityMixingRatio_B::productFunctionSpace
                                              (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("humidity_mixing_ratio_kgkg").functionspace();
}

bool HumidityMixingRatio_B::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering HumidityMixingRatio_B::executeNL function" << std::endl;

    atlas::field::for_each_value(afieldset["humidity_mixing_ratio_kgkg"],
                                 afieldset["humidity_mixing_ratio"],
                                 [&](const double mixr_kgkg, double& mixr) {
        mixr = mixr_kgkg * 1e3;
    });

    oops::Log::trace() << "leaving HumidityMixingRatio_B::executeNL function" << std::endl;

    return true;
}

}  // namespace vader
