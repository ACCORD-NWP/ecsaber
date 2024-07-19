/*
 * (C) Copyright 2024 UCAR
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
#include "vader/recipes/AirPressure.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPressure_B::Name[] = "AirPressure_B";
const oops::JediVariables AirPressure_B::Ingredients{std::vector<std::string>{"base_air_pressure",
                                                             "perturbation_air_pressure"}};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPressure_B> makerAirPressure_B_(AirPressure_B::Name);

// -------------------------------------------------------------------------------------------------

AirPressure_B::AirPressure_B(const AirPressure_BParameters & params,
                             const VaderConfigVars & configVariables) :
                                    configVariables_{configVariables}
{
    oops::Log::trace() << "AirPressure_B::AirPressure_B Starting" << std::endl;
    oops::Log::trace() << "AirPressure_B::AirPressure_B Done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressure_B::name() const {
    return AirPressure_B::Name;
}

// -------------------------------------------------------------------------------------------------

oops::Variable AirPressure_B::product() const {
    return oops::Variable("air_pressure");
}

// -------------------------------------------------------------------------------------------------

oops::JediVariables AirPressure_B::ingredients() const {
    return AirPressure_B::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressure_B::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("perturbation_air_pressure").shape(1);
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressure_B::productFunctionSpace(const atlas::FieldSet & afieldset) const {
    return afieldset.field("perturbation_air_pressure").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressure_B::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "AirPressure_B::executeNL Starting" << std::endl;

    // Get fields
    atlas::Field airPressureBaseF = afieldset.field("base_air_pressure");
    atlas::Field airPressurePertF = afieldset.field("perturbation_air_pressure");
    atlas::Field airPressureF = afieldset.field("air_pressure");

    auto airPressureBase = atlas::array::make_view<double, 2>(airPressureBaseF);
    auto airPressurePert = atlas::array::make_view<double, 2>(airPressurePertF);
    auto airPressure = atlas::array::make_view<double, 2>(airPressureF);

    // Grid dimensions
    size_t h_size = airPressureBaseF.shape(0);
    int v_size = airPressureBaseF.shape(1);

    // Calculate the output variable
    for (int vv = 0; vv < v_size; ++vv) {
      for ( size_t hh = 0; hh < h_size ; ++hh ) {
        airPressure(hh, vv) = airPressureBase(hh, vv) + airPressurePert(hh, vv);
      }
    }

    // Return
    oops::Log::trace() << "AirPressure_B::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
