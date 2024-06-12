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
#include "oops/base/Variable.h"
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
    oops::Variable product() const override {return oops::Variable{"air_pressure_levels"};}
    // oops::Variable product() const override {
    //   return oops::Variable{"air_pressure_extended_up_by_1"};} CCPP Name
    // }

    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;

    bool executeNL(atlas::FieldSet &) override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
