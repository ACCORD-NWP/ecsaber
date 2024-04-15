/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Config.h"

#include "mo/constants.h"
#include "mo/control2analysis_varchange.h"
#include "mo/functions.h"

#include "oops/util/FunctionSpaceHelpers.h"

using atlas::array::make_view;
using atlas::util::Config;
using atlas::idx_t;

namespace mo {

/// \details Calculate the hydrostatic exner pressure (on levels)
///          using air_pressure_minus_one and virtual potential temperature.
void evalHydrostaticExnerLevels(atlas::FieldSet & stateFlds) {
  const auto rpView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto vthetaView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto hexnerView = make_view<double, 2>(stateFlds["hydrostatic_exner_levels"]);

  const idx_t sizeOwned = util::getSizeOwned(stateFlds["hydrostatic_exner_levels"].functionspace());
  const idx_t numLevels = stateFlds["hydrostatic_exner_levels"].shape(1);

  for (idx_t jn = 0; jn < sizeOwned; ++jn) {
    hexnerView(jn, 0) = pow(pView(jn, 0) / constants::p_zero,
      constants::rd_over_cp);
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      hexnerView(jn, jl) = hexnerView(jn, jl-1) -
        (constants::grav * (rpView(jn, jl) - rpView(jn, jl-1))) /
        (constants::cp * vthetaView(jn, jl-1));
    }
  }

  stateFlds["hydrostatic_exner_levels"].set_dirty();
}

/// \details Calculate the hydrostatic pressure (on levels)
///           from hydrostatic exner
void evalHydrostaticPressureLevels(atlas::FieldSet & stateFlds) {
  const auto hexnerView = make_view<double, 2>(stateFlds["hydrostatic_exner_levels"]);
  auto hpView = make_view<double, 2>(stateFlds["hydrostatic_pressure_levels"]);
  const idx_t numLevels = stateFlds["hydrostatic_pressure_levels"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(stateFlds["hydrostatic_pressure_levels"].functionspace());
  for (idx_t jn = 0; jn < sizeOwned; ++jn) {
    for (idx_t jl = 0; jl < numLevels; ++jl) {
       hpView(jn, jl) = constants::p_zero *
         pow(hexnerView(jn, jl), 1.0 / constants::rd_over_cp);
    }
  }
  stateFlds["hydrostatic_pressure_levels"].set_dirty();
}

}  // namespace mo
