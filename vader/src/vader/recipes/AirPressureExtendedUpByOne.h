/*
 * (C) Copyright 2024- UCAR
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

/*! \brief AirPressureExtendedUpByOne_A class defines a recipe for pressure levels as a vertical
           extrapolation of air pressure above model top.
 */
class AirPressureExtendedUpByOne_A : public RecipeBase
{
 public:
    static const char Name[];

    typedef EmptyRecipeParameters Parameters_;

    AirPressureExtendedUpByOne_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override {return Name;}
    std::string product() const override {return "air_pressure_levels";}
    // std::string product() const override {return "air_pressure_extended_up_by_one";} CCPP Name

    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;

    bool executeNL(atlas::FieldSet &) override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
