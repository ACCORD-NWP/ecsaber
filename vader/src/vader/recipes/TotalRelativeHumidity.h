/*
 * (C) Crown Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "vader/RecipeBase.h"

namespace vader
{

// -------------------------------------------------------------------------------------------------

/*! \brief The class 'TotalRelativeHumidity_A' defines a recipe for 'total relative humidity (rht)'
 *
 *  \details This instantiation of RecipeBase produces total relative humidity (rht)
 *           using specific humidity, mass_content_of_cloud_liquid_water_in_atmosphere_layer,
 *           mass_content_of_cloud_ice_in_atmosphere_layer, specific rain (qrain), and
 *           saturated specific humidity (qsat) as inputs.
 *
 */
class TotalRelativeHumidity_A : public RecipeBase
{
 public:
    static const char Name[];

    typedef EmptyRecipeParameters Parameters_;

    TotalRelativeHumidity_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;

    bool executeNL(atlas::FieldSet &) override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
