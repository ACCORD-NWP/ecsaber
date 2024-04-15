
/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "mo/constants.h"
#include "mo/functions.h"
#include "mo/model2geovals_varchange.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;
using atlas::util::Config;

namespace mo {


bool evalParamAParamB(atlas::FieldSet & stateFlds)
{
  oops::Log::trace() << "[evalParamAParamB2()] starting ..." << std::endl;

  std::size_t blindex;
  if (!stateFlds["height"].metadata().has("boundary_layer_index")) {
    oops::Log::error() << "ERROR - data validation failed "
                          "we expect boundary_layer_index value "
                          "in the meta data of the height field" << std::endl;
  }
  stateFlds["height"].metadata().get("boundary_layer_index", blindex);

  const auto heightView = make_view<const double, 2>(stateFlds["height"]);
  const auto heightLevelsView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto pressureLevelsView = make_view<const double, 2>
      (stateFlds["air_pressure_levels_minus_one"]);
  const auto specificHumidityView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  auto param_aView = make_view<double, 2>(stateFlds["surf_param_a"]);
  auto param_bView = make_view<double, 2>(stateFlds["surf_param_b"]);

  // temperature at level above boundary layer
  double t_bl;
  // temperature at model surface height
  double t_msh;

  double exp_pmsh = constants::Lclr * constants::rd / constants::grav;

  const idx_t sizeOwned = util::getSizeOwned(stateFlds["surf_param_a"].functionspace());

  for (idx_t jn = 0; jn < sizeOwned; ++jn) {
    t_bl = (-constants::grav / constants::rd) *
           (heightLevelsView(jn, blindex + 1) - heightLevelsView(jn, blindex)) /
           log(pressureLevelsView(jn, blindex + 1) / pressureLevelsView(jn, blindex));

    t_bl = t_bl / (1.0 + constants::c_virtual * specificHumidityView(jn, blindex));

    t_msh = t_bl + constants::Lclr * (heightView(jn, blindex) - heightLevelsView(jn, 0));

    param_aView(jn, 0) = heightLevelsView(jn, 0) + t_msh / constants::Lclr;
    param_bView(jn, 0) = t_msh / (pow(pressureLevelsView(jn, 0), exp_pmsh) * constants::Lclr);
  }

  stateFlds["surf_param_a"].set_dirty();
  stateFlds["surf_param_b"].set_dirty();

  oops::Log::trace() << "[evalParamAParamB()] ... exit" << std::endl;

  return true;
}

}  // namespace mo
