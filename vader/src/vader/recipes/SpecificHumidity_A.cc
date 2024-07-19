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
#include "vader/recipes/SpecificHumidity.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char SpecificHumidity_A::Name[] = "SpecificHumidity_A";
const oops::JediVariables SpecificHumidity_A::Ingredients{{"humidity_mixing_ratio"}};

// Register the maker
static RecipeMaker<SpecificHumidity_A> makerSpecificHumidity_(SpecificHumidity_A::Name);

SpecificHumidity_A::SpecificHumidity_A(const Parameters_ & params,
                                       const VaderConfigVars & configVariables) :
    configVariables_{configVariables}
{
    oops::Log::trace() << "SpecificHumidity_A::SpecificHumidity_A(params)" << std::endl;
}

std::string SpecificHumidity_A::name() const
{
    return SpecificHumidity_A::Name;
}

oops::Variable SpecificHumidity_A::product() const
{
    return oops::Variable{"humidity_mixing_ratio"};
}

oops::JediVariables SpecificHumidity_A::ingredients() const
{
    return SpecificHumidity_A::Ingredients;
}

size_t SpecificHumidity_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("humidity_mixing_ratio").shape(1);
}

atlas::FunctionSpace SpecificHumidity_A::productFunctionSpace
                                              (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("humidity_mixing_ratio").functionspace();
}

bool SpecificHumidity_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering SpecificHumidity_A::executeNL function" << std::endl;

    // humidity_mixing_ratio in g/kg; specific_humidity in kg/kg
    atlas::field::for_each_value(afieldset["humidity_mixing_ratio"],
                                 afieldset["specific_humidity"],
                                 [&](const double mixr, double& q) {
        q = mixr / (1. + mixr) / 1000.;
    });

    oops::Log::trace() << "leaving SpecificHumidity_A::executeNL function" << std::endl;

    return true;
}

}  // namespace vader
