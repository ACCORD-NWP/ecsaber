/*
 * (C) Crown Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "vader/RecipeBase.h"

namespace vader
{

// -------------------------------------------------------------------------------------------------

/*! \brief The class 'SaturationSpecificHumidity_A' defines a recipe for
 *         'saturation specific humidity (qsat)'
 *
 *  \details This instantiation of RecipeBase produces saturation specific humidity (qsat)
 *           using air pressure, saturation vapour pressure (svp), and air temperature as inputs.
 *
 */
class SaturationSpecificHumidity_A : public RecipeBase
{
 public:
    static const char Name[];

    typedef EmptyRecipeParameters Parameters_;

    SaturationSpecificHumidity_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;

    bool executeNL(atlas::FieldSet &) override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
