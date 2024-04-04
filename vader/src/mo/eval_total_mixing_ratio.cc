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

bool eval_total_mixing_ratio_nl(
        atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] starting ..."
          << std::endl;

  const auto ds_m_v  = make_view<const double, 2>(stateFlds["m_v"]);
  const auto ds_m_ci = make_view<const double, 2>(stateFlds["m_ci"]);
  const auto ds_m_cl = make_view<const double, 2>(stateFlds["m_cl"]);
  const auto ds_m_r  = make_view<const double, 2>(stateFlds["m_r"]);
  auto ds_m_t  = make_view<double, 2>(stateFlds["m_t"]);
  const idx_t sizeOwned =
    util::getSizeOwned(stateFlds["m_t"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < ds_m_t.shape(1); jl++) {
      ds_m_t(jn, jl) = 1.0 + ds_m_v(jn, jl) + ds_m_ci(jn, jl) +
                            ds_m_cl(jn, jl) + ds_m_r(jn, jl);
    }
  }
  stateFlds["m_t"].set_dirty();

  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] ... exit"
          << std::endl;
  return true;
}

// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_tl()] starting ..."
          << std::endl;

  const auto inc_m_v  = make_view<const double, 2>(incFlds["m_v"]);
  const auto inc_m_ci = make_view<const double, 2>(incFlds["m_ci"]);
  const auto inc_m_cl = make_view<const double, 2>(incFlds["m_cl"]);
  const auto inc_m_r  = make_view<const double, 2>(incFlds["m_r"]);
  auto inc_m_t  = make_view<double, 2>(incFlds["m_t"]);
  const idx_t sizeOwned =
    util::getSizeOwned(incFlds["m_t"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < inc_m_t.shape(1); jl++) {
      inc_m_t(jn, jl) = inc_m_v(jn, jl) + inc_m_ci(jn, jl) +
                        inc_m_cl(jn, jl) + inc_m_r(jn, jl);
    }
  }
  incFlds["m_t"].set_dirty();

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

  auto hat_m_v  = make_view<double, 2>(hatFlds["m_v"]);
  auto hat_m_ci = make_view<double, 2>(hatFlds["m_ci"]);
  auto hat_m_cl = make_view<double, 2>(hatFlds["m_cl"]);
  auto hat_m_r  = make_view<double, 2>(hatFlds["m_r"]);
  auto hat_m_t  = make_view<double, 2>(hatFlds["m_t"]);
  const idx_t sizeOwned =
    util::getSizeOwned(hatFlds["m_t"].functionspace());

  atlas_omp_parallel_for(idx_t jn = 0; jn < sizeOwned; jn++) {
    for (idx_t jl = 0; jl < hat_m_t.shape(1); jl++) {
      hat_m_v(jn, jl) += hat_m_t(jn, jl);
      hat_m_ci(jn, jl) += hat_m_t(jn, jl);
      hat_m_cl(jn, jl) += hat_m_t(jn, jl);
      hat_m_r(jn, jl) += hat_m_t(jn, jl);
      hat_m_t(jn, jl) = 0.0;
    }
  }
  hatFlds["m_v"].set_dirty();
  hatFlds["m_ci"].set_dirty();
  hatFlds["m_cl"].set_dirty();
  hatFlds["m_r"].set_dirty();
  hatFlds["m_t"].set_dirty();

  oops::Log::trace()
          << "[eval_total_mixing_ratio_ad()] ... exit"
          << std::endl;
}

// --------------------------------------------------------------------------------------
}  // namespace mo
