/*
 * (C) Copyright 2021-2022 UCAR
 * (C) Crown Copyright 2023-2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_VADER_RECIPES_AIRDENSITYLEVELSMINUSONE_H_
#define SRC_VADER_RECIPES_AIRDENSITYLEVELSMINUSONE_H_

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

// ------------------------------------------------------------------------------------------------

/*! \brief AirDensityLevelsMinusOne_A class defines a recipe for air density
 *
 *  \details This instantiation of RecipeBase produces air density
 *           using virtual_potential_temperature and pressure
 *           (air_pressure_levels_minus_one) as inputs.
 */
class AirDensityLevelsMinusOne_A : public RecipeBase {
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef EmptyRecipeParameters Parameters_;

    AirDensityLevelsMinusOne_A(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

// ------------------------------------------------------------------------------------------------

/*! \brief AirDensityLevelsMinusOne_B class defines a recipe for air density
 *
 *  \details This instantiation of RecipeBase produces air density
 *           using virtual_potential_temperature and pressure (air_pressure_levels)
 *           as inputs.
 */
class AirDensityLevelsMinusOne_B : public RecipeBase {
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef EmptyRecipeParameters Parameters_;

    AirDensityLevelsMinusOne_B(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};


// ------------------------------------------------------------------------------------------------

}  // namespace vader

#endif  // SRC_VADER_RECIPES_AIRDENSITYLEVELSMINUSONE_H_
