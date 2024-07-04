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

class WindReductionFactorAt10m_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(WindReductionFactorAt10m_AParameters, RecipeParametersBase)

 public:
    oops::RequiredParameter<std::string> name{"recipe name", this};
};

/*! \brief WindReductionFactorAt10m_A class defines a recipe to calculate surface wind
 * reduction factor
 *
 *  \details This recipe uses model provided u10/v10 to calculate a ratio against
 *  the winds at the lowest model level for use with the "VertInterp" operator.
 */
class WindReductionFactorAt10m_A : public RecipeBase
{
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef WindReductionFactorAt10m_AParameters Parameters_;

    WindReductionFactorAt10m_A(const Parameters_ &, const VaderConfigVars &);

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
