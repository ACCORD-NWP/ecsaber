/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_surface_wind.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_surface_wind_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_surface_wind_nl()] starting ..." << std::endl;
  const auto uView = make_view<double, 2>(stateFlds["eastward_wind"]);
  const auto vView = make_view<double, 2>(stateFlds["northward_wind"]);
  auto u10mView = make_view<double, 2>(stateFlds["uwind_at_10m"]);
  auto v10mView = make_view<double, 2>(stateFlds["vwind_at_10m"]);
  const idx_t sizeOwned =
    util::getSizeOwned(stateFlds["northward_wind"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; ++jn) {
    u10mView(jn, 0) = uView(jn, 0);
    v10mView(jn, 0) = vView(jn, 0);
  }
  stateFlds["uwind_at_10m"].set_dirty();
  stateFlds["vwind_at_10m"].set_dirty();

  oops::Log::trace() << "[eval_surface_wind_nl()] ... done" << std::endl;
}

void eval_surface_wind_tl(atlas::FieldSet & incFlds) {
  oops::Log::trace() << "[eval_surface_wind_tl()] starting ..." << std::endl;
  const auto uIncView = make_view<double, 2>(incFlds["eastward_wind"]);
  const auto vIncView = make_view<double, 2>(incFlds["northward_wind"]);
  auto u10mIncView = make_view<double, 2>(incFlds["uwind_at_10m"]);
  auto v10mIncView = make_view<double, 2>(incFlds["vwind_at_10m"]);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["uwind_at_10m"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; ++jn) {
      u10mIncView(jn, 0) = uIncView(jn, 0);
      v10mIncView(jn, 0) = vIncView(jn, 0);
  }
  incFlds["uwind_at_10m"].set_dirty();
  incFlds["vwind_at_10m"].set_dirty();

  oops::Log::trace() << "[eval_surface_wind_tl()] ... done" << std::endl;
}

void eval_surface_wind_ad(atlas::FieldSet & hatFlds) {
  oops::Log::trace() << "[eval_surface_wind_ad()] starting ..." << std::endl;
  auto uHatView = make_view<double, 2>(hatFlds["eastward_wind"]);
  auto vHatView = make_view<double, 2>(hatFlds["northward_wind"]);
  auto u10mHatView = make_view<double, 2>(hatFlds["uwind_at_10m"]);
  auto v10mHatView = make_view<double, 2>(hatFlds["vwind_at_10m"]);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["uwind_at_10m"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; ++jn) {
      uHatView(jn, 0) += u10mHatView(jn, 0);
      vHatView(jn, 0) += v10mHatView(jn, 0);
      u10mHatView(jn, 0) = 0.0;
      v10mHatView(jn, 0) = 0.0;
  }
  hatFlds["eastward_wind"].set_dirty();
  hatFlds["northward_wind"].set_dirty();
  hatFlds["uwind_at_10m"].set_dirty();
  hatFlds["vwind_at_10m"].set_dirty();

  oops::Log::trace() << "[eval_surface_wind_ad()] ... done" << std::endl;
}

}  // namespace mo
