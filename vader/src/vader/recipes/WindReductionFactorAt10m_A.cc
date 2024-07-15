/*
 * (C) Copyright 2024  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "oops/util/Logger.h"
#include "vader/recipes/WindReductionFactorAt10m.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char WindReductionFactorAt10m_A::Name[] = "WindReductionFactorAt10m_A";
const oops::JediVariables WindReductionFactorAt10m_A::Ingredients{
      std::vector<std::string>{"surface_eastward_wind", "surface_northward_wind",
                               "eastward_wind", "northward_wind"}};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<WindReductionFactorAt10m_A>
       makerWindReductionFactorAt10m_A_(WindReductionFactorAt10m_A::Name);

// -------------------------------------------------------------------------------------------------

WindReductionFactorAt10m_A::WindReductionFactorAt10m_A(
                            const WindReductionFactorAt10m_AParameters & params,
                            const VaderConfigVars & configVariables):
                            configVariables_{configVariables}
{
    oops::Log::trace() << "WindReductionFactorAt10m_A::WindReductionFactorAt10m_A Starting"
                       << std::endl;
    oops::Log::trace() << "WindReductionFactorAt10m_A::WindReductionFactorAt10m_A Done"
                       << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string WindReductionFactorAt10m_A::name() const {
    return WindReductionFactorAt10m_A::Name;
}

// -------------------------------------------------------------------------------------------------

    // variable name used in yaml file
oops::Variable WindReductionFactorAt10m_A::product() const {
    return oops::Variable("wind_reduction_factor_at_10m");
}

// -------------------------------------------------------------------------------------------------

oops::JediVariables WindReductionFactorAt10m_A::ingredients() const {
    return WindReductionFactorAt10m_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t WindReductionFactorAt10m_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("surface_eastward_wind").shape(1);
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace WindReductionFactorAt10m_A::productFunctionSpace(const atlas::FieldSet
& afieldset) const {
    return afieldset.field("surface_eastward_wind").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool WindReductionFactorAt10m_A::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "WindReductionFactorAt10m_A::executeNL Starting" << std::endl;

    // Get the fields
    atlas::Field uu = afieldset.field("eastward_wind");
    atlas::Field vv = afieldset.field("northward_wind");
    atlas::Field u10m = afieldset.field("surface_eastward_wind");
    atlas::Field v10m = afieldset.field("surface_northward_wind");
    // output:
    atlas::Field f10m = afieldset.field("wind_reduction_factor_at_10m");

    // Set the units
    f10m.metadata().set("units", "none");

    // Set the array views to manipulate the data
    auto uu_view = atlas::array::make_view<double, 2>(uu);
    auto vv_view = atlas::array::make_view<double, 2>(vv);
    auto u10m_view = atlas::array::make_view<double, 2>(u10m);
    auto v10m_view = atlas::array::make_view<double, 2>(v10m);
    auto f10m_view = atlas::array::make_view<double, 2>(f10m);

    // Get the grid size
    const int gridSize = uu.shape(0);    // geom%nCells
    const int nLevel = uu.shape(1) - 1;  // Reduce by 1 since index begins at 0

    // for atlas fields, vertical is: top --> bottom.
    //
    for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
        f10m_view(jNode, 0) = sqrt(pow(u10m_view(jNode, 0), 2) + pow(v10m_view(jNode, 0), 2));
        if (f10m_view(jNode, 0) > 0) {
          f10m_view(jNode, 0) = f10m_view(jNode, 0)/sqrt(pow(uu_view(jNode, nLevel), 2)
                               + pow(vv_view(jNode, nLevel), 2));
        } else {
          f10m_view(jNode, 0) = 1.0;
        }
    }

    // Return
    oops::Log::trace() << "WindReductionFactorAt10m_A::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
