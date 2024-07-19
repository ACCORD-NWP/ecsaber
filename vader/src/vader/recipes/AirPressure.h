/*
 * (C) Copyright 2024 UCAR
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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader
{

// -------------------------------------------------------------------------------------------------

class AirPressure_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirPressure_AParameters, RecipeParametersBase)

 public:
    oops::RequiredParameter<std::string> name{"recipe name", this};
};

class AirPressure_BParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirPressure_BParameters, RecipeParametersBase)

 public:
    oops::RequiredParameter<std::string> name{"recipe name", this};
};
/*! \brief AirPressure_A class defines a recipe for pressure levels from pressure
           thickness.
 *
 *  \details This recipe uses pressure at the interfaces, along with the Phillips method to
 *           compute pressure at the mid points. It does not provide TL/AD algorithms.
 */
class AirPressure_A : public RecipeBase
{
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef AirPressure_AParameters Parameters_;

    AirPressure_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

/*! \brief AirPressure_B class defines a recipe for pressure levels from perturbation air pressure
           and base pressure (for WRF).
 *
 *  \details This recipe calculates AirPressure at mid points by summation of perturbation air pressure
 *           and base air pressure. It does not provide TL/AD algorithms.
 */
class AirPressure_B : public RecipeBase
{
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef AirPressure_BParameters Parameters_;

    AirPressure_B(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
