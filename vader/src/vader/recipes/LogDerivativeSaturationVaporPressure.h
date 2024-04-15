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

/*! \brief The class 'LogDerivativeSaturationVaporPressure_A' defines a recipe for
 *         'derivative of ln water vapor partial pressure assuming saturation at interface
 *         wrt air temperature at interface (dlsvpdT)'
 *
 *  \details This instantiation of RecipeBase produces derivative of ln water vapor partial
 *           pressure assuming saturation at interface wrt air temperature at interface (dlsvpdT)
 *           using air temperature as input.
 *
 */
class LogDerivativeSaturationVaporPressure_A : public RecipeBase
{
 public:
    static const char Name[];

    typedef EmptyRecipeParameters Parameters_;

    LogDerivativeSaturationVaporPressure_A(const Parameters_ &,
                                           const VaderConfigVars &);

    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;

    bool executeNL(atlas::FieldSet &) override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
