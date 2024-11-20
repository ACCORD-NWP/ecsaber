/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_ratio.h"
#include "mo/eval_total_mixing_ratio.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

namespace mo {

using atlas::array::make_view;
using atlas::idx_t;

// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_nl(
        atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] starting ..."
          << std::endl;

  const auto ds_m_v  = make_view<const double, 2>(stateFlds[
                                                   "water_vapor_mixing_ratio_wrt_dry_air"]);
  const auto ds_m_ci = make_view<const double, 2>(stateFlds["cloud_ice_mixing_ratio_wrt_dry_air"]);
  const auto ds_m_cl = make_view<const double, 2>(stateFlds[
                                                   "cloud_liquid_water_mixing_ratio_wrt_dry_air"]);
  const auto ds_m_r  = make_view<const double, 2>(stateFlds["rain_mixing_ratio_wrt_dry_air"]);
  auto ds_m_t  = make_view<double, 2>(stateFlds["total_water_mixing_ratio_wrt_dry_air"]);
  const idx_t sizeOwned =
    util::getSizeOwned(stateFlds["total_water_mixing_ratio_wrt_dry_air"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < ds_m_t.shape(1); jl++) {
      ds_m_t(jn, jl) = 1.0 + ds_m_v(jn, jl) + ds_m_ci(jn, jl) +
                            ds_m_cl(jn, jl) + ds_m_r(jn, jl);
    }
  }
  stateFlds["total_water_mixing_ratio_wrt_dry_air"].set_dirty();

  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] ... exit"
          << std::endl;
}

// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_tl()] starting ..."
          << std::endl;

  const auto inc_m_v  = make_view<const double, 2>(incFlds["water_vapor_mixing_ratio_wrt_dry_air"]);
  const auto inc_m_ci = make_view<const double, 2>(incFlds["cloud_ice_mixing_ratio_wrt_dry_air"]);
  const auto inc_m_cl = make_view<const double, 2>(incFlds[
                                                    "cloud_liquid_water_mixing_ratio_wrt_dry_air"]);
  const auto inc_m_r  = make_view<const double, 2>(incFlds["rain_mixing_ratio_wrt_dry_air"]);
  auto inc_m_t  = make_view<double, 2>(incFlds["total_water_mixing_ratio_wrt_dry_air"]);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["total_water_mixing_ratio_wrt_dry_air"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < inc_m_t.shape(1); jl++) {
      inc_m_t(jn, jl) = inc_m_v(jn, jl) + inc_m_ci(jn, jl) +
                        inc_m_cl(jn, jl) + inc_m_r(jn, jl);
    }
  }
  incFlds["total_water_mixing_ratio_wrt_dry_air"].set_dirty();

  oops::Log::trace()
          << "[eval_total_mixing_ratio_tl()] ... exit"
          << std::endl;
}
// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_ad()] starting ..."
          << std::endl;

  auto hat_m_v  = make_view<double, 2>(hatFlds["water_vapor_mixing_ratio_wrt_dry_air"]);
  auto hat_m_ci = make_view<double, 2>(hatFlds["cloud_ice_mixing_ratio_wrt_dry_air"]);
  auto hat_m_cl = make_view<double, 2>(hatFlds["cloud_liquid_water_mixing_ratio_wrt_dry_air"]);
  auto hat_m_r  = make_view<double, 2>(hatFlds["rain_mixing_ratio_wrt_dry_air"]);
  auto hat_m_t  = make_view<double, 2>(hatFlds["total_water_mixing_ratio_wrt_dry_air"]);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["total_water_mixing_ratio_wrt_dry_air"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < hat_m_t.shape(1); jl++) {
      hat_m_v(jn, jl) += hat_m_t(jn, jl);
      hat_m_ci(jn, jl) += hat_m_t(jn, jl);
      hat_m_cl(jn, jl) += hat_m_t(jn, jl);
      hat_m_r(jn, jl) += hat_m_t(jn, jl);
      hat_m_t(jn, jl) = 0.0;
    }
  }
  hatFlds["water_vapor_mixing_ratio_wrt_dry_air"].set_dirty();
  hatFlds["cloud_ice_mixing_ratio_wrt_dry_air"].set_dirty();
  hatFlds["cloud_liquid_water_mixing_ratio_wrt_dry_air"].set_dirty();
  hatFlds["rain_mixing_ratio_wrt_dry_air"].set_dirty();
  hatFlds["total_water_mixing_ratio_wrt_dry_air"].set_dirty();

  oops::Log::trace()
          << "[eval_total_mixing_ratio_ad()] ... exit"
          << std::endl;
}

// --------------------------------------------------------------------------------------
}  // namespace mo
