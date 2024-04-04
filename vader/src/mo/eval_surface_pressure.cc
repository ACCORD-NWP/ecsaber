/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_surface_pressure.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_surface_pressure_tl(atlas::FieldSet & incFlds) {
  oops::Log::trace() << "[eval_surface_pressure_tl()] starting ..." << std::endl;
  const auto pIncView = make_view<double, 2>(incFlds["air_pressure_levels"]);
  auto pstarIncView = make_view<double, 2>(incFlds["surface_pressure"]);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["surface_pressure"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; ++jn) {
      pstarIncView(jn, 0) = pIncView(jn, 0);
  }
  incFlds["surface_pressure"].set_dirty();
  oops::Log::trace() << "[eval_surface_pressure_tl()] ... done" << std::endl;
}

void eval_surface_pressure_ad(atlas::FieldSet & hatFlds) {
  oops::Log::trace() << "[eval_surface_pressure_ad()] starting ..." << std::endl;
  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto pstarHatView = make_view<double, 2>(hatFlds["surface_pressure"]);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["surface_pressure"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; ++jn) {
    pHatView(jn, 0) += pstarHatView(jn, 0);
    pstarHatView(jn, 0) = 0.0;
  }
  hatFlds["air_pressure_levels"].set_dirty();
  hatFlds["surface_pressure"].set_dirty();
  oops::Log::trace() << "[eval_surface_pressure_ad()] ... done" << std::endl;
}

}  // namespace mo
