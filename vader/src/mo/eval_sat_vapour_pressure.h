/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field/FieldSet.h"


namespace mo {


/// \brief function to evaluate the several types of 'saturation water pressure'
///
///  input: Atlas data structure (FieldSet) containing 'air temperature'
/// output: 'svp' or 'dlsvpdT'
///
void eval_sat_vapour_pressure_variant_nl(const std::string field_out,
                                         atlas::FieldSet & fields);

/// \brief function to evaluate the 'saturation water pressure (svp)'
void eval_sat_vapour_pressure_nl(atlas::FieldSet & fields);

/// \brief function to evaluate the 'derivative of ln water vapor partial pressure
/// assuming saturation at interface wrt air temperature at interface (dlsvpdT)'
void eval_derivative_ln_svp_wrt_temperature_nl(atlas::FieldSet & fields);


}  // namespace mo
