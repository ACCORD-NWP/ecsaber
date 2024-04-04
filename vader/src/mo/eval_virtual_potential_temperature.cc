/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/constants.h"
#include "mo/eval_virtual_potential_temperature.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_nl()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  auto vThetaView = make_view<double, 2>(stateFlds["virtual_potential_temperature"]);

  const idx_t n_levels = stateFlds["virtual_potential_temperature"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(stateFlds["virtual_potential_temperature"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < n_levels; jl++) {
      vThetaView(jn, jl) = thetaView(jn, jl)
        * (1.0 + constants::c_virtual * qView(jn, jl));
    }
  }
  stateFlds["virtual_potential_temperature"].set_dirty();

  oops::Log::trace() << "[eval_virtual_potential_temperature_nl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_tl(atlas::FieldSet & incFlds,
                                           const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_tl()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  const auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto vThetaIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const idx_t n_levels = incFlds["virtual_potential_temperature"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["virtual_potential_temperature"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < n_levels; jl++) {
      vThetaIncView(jn, jl) = thetaView(jn, jl) * constants::c_virtual * qIncView(jn, jl) +
                              thetaIncView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl));
    }
  }
  incFlds["virtual_potential_temperature"].set_dirty();
  oops::Log::trace() << "[eval_virtual_potential_temperature_tl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_ad(atlas::FieldSet & hatFlds,
                                           const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_ad()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto vThetaHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);

  const idx_t n_levels = hatFlds["virtual_potential_temperature"].shape(1);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["virtual_potential_temperature"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < n_levels; jl++) {
      qHatView(jn, jl) += thetaView(jn, jl) * constants::c_virtual * vThetaHatView(jn, jl);
      thetaHatView(jn, jl) +=  vThetaHatView(jn, jl)
            * (1.0 + constants::c_virtual * qView(jn, jl));
      vThetaHatView(jn, jl) = 0;
    }
  }
  hatFlds["specific_humidity"].set_dirty();
  hatFlds["potential_temperature"].set_dirty();
  hatFlds["virtual_potential_temperature"].set_dirty();

  oops::Log::trace() << "[eval_virtual_potential_temperature_ad()] ... done" << std::endl;
}

}  // namespace mo
