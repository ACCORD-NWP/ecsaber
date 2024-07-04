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
#include "vader/recipes/AirPressureAtInterface.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPressureAtInterface_C::Name[] = "AirPressureAtInterface_C";
const oops::JediVariables AirPressureAtInterface_C::Ingredients{
      std::vector<std::string>{"surface_pressure", "air_pressure",
                               "geometric_height_of_layer_interfaces"}};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPressureAtInterface_C>
       makerAirPressureAtInterface_C_(AirPressureAtInterface_C::Name);

// -------------------------------------------------------------------------------------------------

AirPressureAtInterface_C::AirPressureAtInterface_C(
                          const AirPressureAtInterface_CParameters & params,
                          const VaderConfigVars & configVariables):
                          configVariables_{configVariables}
{
    oops::Log::trace() << "AirPressureAtInterface_C::AirPressureAtInterface_C Starting"
                       << std::endl;
    oops::Log::trace() << "AirPressureAtInterface_C::AirPressureAtInterface_C Done"
                       << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureAtInterface_C::name() const {
    return AirPressureAtInterface_C::Name;
}

// -------------------------------------------------------------------------------------------------

oops::Variable AirPressureAtInterface_C::product() const {
    return oops::Variable("air_pressure_levels");
}

// -------------------------------------------------------------------------------------------------

oops::JediVariables AirPressureAtInterface_C::ingredients() const {
    return AirPressureAtInterface_C::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressureAtInterface_C::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("geometric_height_of_layer_interfaces").shape(1);
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressureAtInterface_C::productFunctionSpace(const atlas::FieldSet
& afieldset) const {
    return afieldset.field("geometric_height_of_layer_interfaces").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressureAtInterface_C::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "AirPressureAtInterface_C::executeNL Starting" << std::endl;

    // Get the input fields
    atlas::Field ps = afieldset.field("surface_pressure");
    atlas::Field p = afieldset.field("air_pressure");
    atlas::Field zgrid = afieldset.field("geometric_height_of_layer_interfaces");
    // output:
    atlas::Field prsi = afieldset.field("air_pressure_levels");

    // Get and set the units
    std::string p_units, prsi_units;
    p.metadata().get("units", p_units);
    prsi.metadata().set("units", p_units);

    // Set the array views to manipulate the data
    auto zgrid_view = atlas::array::make_view<double, 2>(zgrid);
    auto p_view = atlas::array::make_view<double, 2>(p);
    auto ps_view = atlas::array::make_view<double, 2>(ps);
    auto prsi_view = atlas::array::make_view<double, 2>(prsi);

    // Get the grid size
    const int gridSize = zgrid.shape(0);    // geom%nCells
    const int nLevel = zgrid.shape(1) - 1;  // Reduce by 1 since index begins at 0

    // for atlas, vertical is top --> bottom. Start from bottom to top
    // Set pressure at the bottom to surface pressure
    for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
        prsi_view(jNode, nLevel) = ps_view(jNode, 0);
    }

    // Interpolate linearly pressure log(p) to geometrical height levels (zgrid)
    // (follow MPAS-JEDI subroutine pressure_half_to_full). start from bottom to top
    //
    for (int level = nLevel-1; level >=1; --level) {
        for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
            double w1 = (zgrid_view(jNode, level) - zgrid_view(jNode, level+1))
                        /(zgrid_view(jNode, level-1) - zgrid_view(jNode, level+1));
            double w2 = 1.0 - w1;
            prsi_view(jNode, level) = exp(w1*log(p_view(jNode, level))
                        + w2*log(p_view(jNode, level+1)));
        }
            std::cout << std::scientific << prsi_view(2, level) << '\n';
    }

    // Extrapolate for the top level
    int level = 0;
    for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
        double z0 = zgrid_view(jNode, level);
        double z1 = 0.5*(zgrid_view(jNode, level) + zgrid_view(jNode, level+1));
        double z2 = 0.5*(zgrid_view(jNode, level+1) + zgrid_view(jNode, level+2));
        double w1 = (z0-z2)/(z1-z2);
        double w2 = 1.0 - w1;
        prsi_view(jNode, level) = exp(w1*log(p_view(jNode, level+1))
                                  + w2*log(p_view(jNode, level+2)) );
        std::cout << std::scientific << prsi_view(jNode, level) << '\n';
    }

    // Return
    oops::Log::trace() << "AirPressureAtInterface_C::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
