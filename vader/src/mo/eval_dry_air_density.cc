/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"

#include "mo/constants.h"
#include "mo/eval_dry_air_density.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

// -------------------------------------------------------------------------------------------------

/// \details Calculate the dry air density from virtual potential temperature
///          and air pressure (using air_pressure_levels_minus_one)
void eval_dry_air_density_from_pressure_levels_minus_one_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = stateFlds["dry_air_density_levels_minus_one"].shape(1);

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    rhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptView(jn, 0));

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoView(jn, jl) = pView(jn, jl) *
        pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp)
        * h_minus_hm1 / (constants::rd * (
        h_minus_hl * vptView(jn, jl-1) + hl_minus_hm1 * vptView(jn, jl)));
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                      const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels_minus_one"]);
  auto vptIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto rhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = incFlds["dry_air_density_levels_minus_one"].shape(1);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoIncView.shape(0); ++jn) {
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoIncView(jn, jl) = rhoView(jn, jl) * ( (1.0 - constants::rd_over_cp) *
        pIncView(jn, jl) / pView(jn, jl) -
        ((hl_minus_hm1 * vptIncView(jn, jl) + h_minus_hl * vptIncView(jn, jl-1)) /
         (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1))));
    }

    rhoIncView(jn, 0) = rhoView(jn, 0) * (
      (1.0 - constants::rd_over_cp) * pIncView(jn, 0) / pView(jn, 0) -
      vptIncView(jn, 0)/ vptView(jn, 0));
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                                     const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_ad()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels_minus_one"]);
  auto vptHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto rhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoHatView.shape(0); ++jn) {
    pHatView(jn, 0) += rhoView(jn, 0) * rhoHatView(jn, 0) * (1.0 - constants::rd_over_cp) /
      pView(jn, 0);
    vptHatView(jn, 0) -= rhoView(jn, 0) * rhoHatView(jn, 0) /
      vptView(jn, 0);
    rhoHatView(jn, 0) = 0.0;

    for (idx_t jl = hatFlds["dry_air_density_levels_minus_one"].shape(1)-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += rhoView(jn, jl) * rhoHatView(jn, jl) * (1.0 - constants::rd_over_cp) /
        pView(jn, jl);
      vptHatView(jn, jl) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        hl_minus_hm1 /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      vptHatView(jn, jl-1) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        h_minus_hl /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      rhoHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_ad()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

/// \details Calculate the dry air density from virtual potential temperature
///          and air pressure (using air_pressure_levels)
void eval_dry_air_density_from_pressure_levels_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = stateFlds["dry_air_density_levels_minus_one"].shape(1);

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    rhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptView(jn, 0));

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoView(jn, jl) = pView(jn, jl) *
        pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp)
        * h_minus_hm1 / (constants::rd * (
        h_minus_hl * vptView(jn, jl-1) + hl_minus_hm1 * vptView(jn, jl)));
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_tl()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels"]);
  auto vptIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto rhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = incFlds["dry_air_density_levels_minus_one"].shape(1);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoIncView.shape(0); ++jn) {
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoIncView(jn, jl) = rhoView(jn, jl) * ( (1.0 - constants::rd_over_cp) *
        pIncView(jn, jl) / pView(jn, jl) -
        ((hl_minus_hm1 * vptIncView(jn, jl) + h_minus_hl * vptIncView(jn, jl-1)) /
         (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1))));
    }

    rhoIncView(jn, 0) = rhoView(jn, 0) * (
      (1.0 - constants::rd_over_cp) * pIncView(jn, 0) / pView(jn, 0) -
      vptIncView(jn, 0)/ vptView(jn, 0));
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_ad()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto vptHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto rhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoHatView.shape(0); ++jn) {
    pHatView(jn, 0) += rhoView(jn, 0) * rhoHatView(jn, 0) * (1.0 - constants::rd_over_cp) /
      pView(jn, 0);
    vptHatView(jn, 0) -= rhoView(jn, 0) * rhoHatView(jn, 0) /
      vptView(jn, 0);
    rhoHatView(jn, 0) = 0.0;

    for (idx_t jl = hatFlds["dry_air_density_levels_minus_one"].shape(1)-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += rhoView(jn, jl) * rhoHatView(jn, jl) * (1.0 - constants::rd_over_cp) /
        pView(jn, jl);
      vptHatView(jn, jl) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        hl_minus_hm1 /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      vptHatView(jn, jl-1) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        h_minus_hl /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      rhoHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_ad()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

/// \details Calculate the correct dry air density from potential temperature,
///          specific humidity, mass_content_of_cloud_liquid_water_in_atmosphere_layer,
///          mass_content_of_cloud_ice_in_atmosphere_layer
///          and air pressure (using air_pressure_levels)
void eval_new_dry_air_density_from_pressure_levels_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto qclView = make_view<const double, 2>(
                         stateFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qcfView = make_view<const double, 2>(
                         stateFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto dryrhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydens_jlm1;
  double vptdrydens_intp;

  idx_t nlevels = stateFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    // vptdrydens is the virtual potential temperature for dry density (theta_vdd)
    // so that the dry air density (rho_dd) can be written as rho_dd = p / (Rd Pi theta_vdd)
    vptdrydens = ptView(jn, 0) * (1.0 + constants::c_virtual * qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0)) / (1.0 - qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0));
    dryrhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptdrydens);

    for (idx_t jl = 1; jl < nlevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens_jlm1 = vptdrydens;
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl)) / (1.0 - qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydens_intp = (h_minus_hl * vptdrydens_jlm1 + hl_minus_hm1 * vptdrydens) / h_minus_hm1;
      dryrhoView(jn, jl) = pView(jn, jl) *
                           pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp) /
                           (constants::rd * vptdrydens_intp);
    }
  }
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_new_dry_air_density_from_pressure_levels_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_tl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto qclView = make_view<const double, 2>(
                         stateFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qcfView = make_view<const double, 2>(
                         stateFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  const auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  const auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels"]);
  const auto ptIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  const auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  const auto qclIncView = make_view<const double, 2>(
                            incFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qcfIncView = make_view<const double, 2>(
                            incFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto dryrhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydensInc;
  double vptdrydens_jlm1;
  double vptdrydensInc_jlm1;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensInc_intp_times_h_minus_hm1;

  idx_t nlevels = incFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < dryrhoIncView.shape(0); ++jn) {
    vptdrydens = ptView(jn, 0) * (1.0 + constants::c_virtual * qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0)) / (1.0 - qView(jn, 0)
                 - qclView(jn, 0) - qcfView(jn, 0));
    vptdrydensInc = ((1.0 + constants::c_virtual * qView(jn, 0)
                      - qclView(jn, 0) - qcfView(jn, 0)) * ptIncView(jn, 0)
                     + ((1.0 + constants::c_virtual) * qIncView(jn, 0)
                        * (1.0 - qclView(jn, 0) - qcfView(jn, 0)) * ptView(jn, 0)
                        / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)))
                     + (qclIncView(jn, 0) *  qView(jn, 0)
                       * (1.0 + constants::c_virtual) * ptView(jn, 0))
                       / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                     + (qcfIncView(jn, 0) *  qView(jn, 0)
                       * (1.0 + constants::c_virtual) * ptView(jn, 0))
                       / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)))
                    / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0));
    dryrhoIncView(jn, 0) = dryrhoView(jn, 0) * ((1.0 - constants::rd_over_cp)
                           * pIncView(jn, 0) / pView(jn, 0)
                           - vptdrydensInc / vptdrydens);

    for (idx_t jl = 1; jl < nlevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens_jlm1 = vptdrydens;
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl)) / (1.0 - qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydensInc_jlm1 = vptdrydensInc;
      vptdrydensInc = ((1.0 + constants::c_virtual * qView(jn, jl)
                        - qclView(jn, jl) - qcfView(jn, jl)) * ptIncView(jn, jl)
                       + ((1.0 + constants::c_virtual) * qIncView(jn, jl)
                          * (1.0 - qclView(jn, jl) - qcfView(jn, jl)) * ptView(jn, jl)
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                       + (qclIncView(jn, jl) *  qView(jn, jl)
                          * (1.0 + constants::c_virtual) * ptView(jn, jl))
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                       + (qcfIncView(jn, jl) *  qView(jn, jl)
                          * (1.0 + constants::c_virtual) * ptView(jn, jl))
                          / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                      / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensInc_intp_times_h_minus_hm1 = hl_minus_hm1 * vptdrydensInc +
                                              h_minus_hl * vptdrydensInc_jlm1;
      dryrhoIncView(jn, jl) = dryrhoView(jn, jl) * ( (1.0 - constants::rd_over_cp)
                              * pIncView(jn, jl) / pView(jn, jl)
                              - vptdrydensInc_intp_times_h_minus_hm1
                              / vptdrydens_intp_times_h_minus_hm1);
    }
  }
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_new_dry_air_density_from_pressure_levels_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_ad()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto qclView = make_view<const double, 2>(
                 stateFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qcfView = make_view<const double, 2>(
                       stateFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  const auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto ptHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto qclHatView = make_view<double, 2>(
                    hatFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfHatView = make_view<double, 2>(hatFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto dryrhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensHat;
  double vptdrydens_jlm1;
  double vptdrydensHat_jlm1;

  idx_t nlevels = hatFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < dryrhoHatView.shape(0); ++jn) {
    for (idx_t jl = nlevels-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += dryrhoView(jn, jl) * (1.0 - constants::rd_over_cp) /
                          pView(jn, jl) * dryrhoHatView(jn, jl);
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl)) / (1.0 - qView(jn, jl)
                   - qclView(jn, jl) - qcfView(jn, jl));
      vptdrydens_jlm1 = ptView(jn, jl-1) * (1.0 + constants::c_virtual * qView(jn, jl-1)
                        - qclView(jn, jl-1) - qcfView(jn, jl-1)) / (1.0 - qView(jn, jl-1)
                        - qclView(jn, jl-1) - qcfView(jn, jl-1));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensHat = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * hl_minus_hm1
                      / vptdrydens_intp_times_h_minus_hm1;
      vptdrydensHat_jlm1 = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * h_minus_hl /
                           vptdrydens_intp_times_h_minus_hm1;
      ptHatView(jn, jl) += (1.0 + constants::c_virtual * qView(jn, jl)
                            - qclView(jn, jl) - qcfView(jn, jl))
                           / (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                           * vptdrydensHat;
      qHatView(jn, jl) += (1.0 + constants::c_virtual) * ptView(jn, jl)
                           * (1.0 - qclView(jn, jl) - qcfView(jn, jl))
                           / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                               * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                           * vptdrydensHat;
      qclHatView(jn, jl) += ptView(jn, jl) *  qView(jn, jl) * (1.0 + constants::c_virtual)
                            / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                                * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                            * vptdrydensHat;
      qcfHatView(jn, jl) += ptView(jn, jl) *  qView(jn, jl) * (1.0 + constants::c_virtual)
                            / ((1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl))
                                * (1.0 - qView(jn, jl) - qclView(jn, jl) - qcfView(jn, jl)))
                            * vptdrydensHat;
      ptHatView(jn, jl-1) += (1.0 + constants::c_virtual * qView(jn, jl-1)
                              - qclView(jn, jl-1) - qcfView(jn, jl-1))
                             / (1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                             * vptdrydensHat_jlm1;
      qHatView(jn, jl-1) += (1.0 + constants::c_virtual) * ptView(jn, jl-1)
                            * (1.0 - qclView(jn, jl-1) - qcfView(jn, jl-1))
                            / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1)))
                            * vptdrydensHat_jlm1;
      qclHatView(jn, jl-1) += ptView(jn, jl-1) *  qView(jn, jl-1) * (1.0 + constants::c_virtual)
                              / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                  * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1)
                                     - qcfView(jn, jl-1)))
                              * vptdrydensHat_jlm1;
      qcfHatView(jn, jl-1) += ptView(jn, jl-1) *  qView(jn, jl-1) * (1.0 + constants::c_virtual)
                              / ((1.0 - qView(jn, jl-1) - qclView(jn, jl-1) - qcfView(jn, jl-1))
                                  * (1.0 - qView(jn, jl-1) - qclView(jn, jl-1)
                                     - qcfView(jn, jl-1)))
                              * vptdrydensHat_jlm1;
      dryrhoHatView(jn, jl) = 0.0;
    }
    pHatView(jn, 0) += dryrhoView(jn, 0) * ((1.0 - constants::rd_over_cp) /
                       pView(jn, 0)) * dryrhoHatView(jn, 0);
    vptdrydens = ptView(jn, 0) *
                 (1.0 + constants::c_virtual * qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)) /
                 (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0));
    vptdrydensHat = - dryrhoView(jn, 0) * dryrhoHatView(jn, 0) / vptdrydens;
    ptHatView(jn, 0) += (1.0 + constants::c_virtual * qView(jn, 0)
                         - qclView(jn, 0) - qcfView(jn, 0))
                        / (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0)) * vptdrydensHat;
    qHatView(jn, 0) += (1.0 + constants::c_virtual) * ptView(jn, 0)
                       * (1.0 - qclView(jn, 0) - qcfView(jn, 0))
                       / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                       * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    qclHatView(jn, 0) += ptView(jn, 0) *  qView(jn, 0) * (1.0 + constants::c_virtual)
                         / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                         * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    qcfHatView(jn, 0) += ptView(jn, 0) *  qView(jn, 0) * (1.0 + constants::c_virtual)
                         / ((1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))
                         * (1.0 - qView(jn, 0) - qclView(jn, 0) - qcfView(jn, 0))) * vptdrydensHat;
    dryrhoHatView(jn, 0) = 0.0;
  }
  oops::Log::trace() << "[eval_new_dry_air_density_from_pressure_levels_ad()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

/// \details Calculate the correct dry air density from potential temperature,
///          specific humidity and air pressure (using air_pressure_levels)
void eval_new_dry_air_density_without_condensate_from_pressure_levels_nl(atlas::FieldSet
                                                                         & stateFlds) {
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto dryrhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydens_jlm1;
  double vptdrydens_intp;

  idx_t nlevels = stateFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    // vptdrydens is the virtual potential temperature for dry density (theta_vdd)
    // so that the dry air density (rho_dd) can be written as rho_dd = p / (Rd Pi theta_vdd)
    vptdrydens = ptView(jn, 0) * (1.0 + constants::c_virtual * qView(jn, 0))
                 / (1.0 - qView(jn, 0));
    dryrhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptdrydens);

    for (idx_t jl = 1; jl < nlevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens_jlm1 = vptdrydens;
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl))
                   / (1.0 - qView(jn, jl));
      vptdrydens_intp = (h_minus_hl * vptdrydens_jlm1 + hl_minus_hm1 * vptdrydens) / h_minus_hm1;
      dryrhoView(jn, jl) = pView(jn, jl) *
                           pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp) /
                           (constants::rd * vptdrydens_intp);
    }
  }
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_new_dry_air_density_without_condensate_from_pressure_levels_tl(atlas::FieldSet & incFlds,
                                                                         const atlas::FieldSet
                                                                         & stateFlds) {
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_tl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  const auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  const auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels"]);
  const auto ptIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  const auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  auto dryrhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydensInc;
  double vptdrydens_jlm1;
  double vptdrydensInc_jlm1;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensInc_intp_times_h_minus_hm1;

  idx_t nlevels = incFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < dryrhoIncView.shape(0); ++jn) {
    vptdrydens = ptView(jn, 0) * (1.0 + constants::c_virtual * qView(jn, 0))
                 / (1.0 - qView(jn, 0));
    vptdrydensInc = ((1.0 + constants::c_virtual * qView(jn, 0)) * ptIncView(jn, 0)
                     + ((1.0 + constants::c_virtual) * qIncView(jn, 0)
                     * ptView(jn, 0)
                     / (1.0 - qView(jn, 0))))
                    / (1.0 - qView(jn, 0));
    dryrhoIncView(jn, 0) = dryrhoView(jn, 0) * ((1.0 - constants::rd_over_cp)
                           * pIncView(jn, 0) / pView(jn, 0)
                           - vptdrydensInc / vptdrydens);

    for (idx_t jl = 1; jl < nlevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      vptdrydens_jlm1 = vptdrydens;
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl))
                   / (1.0 - qView(jn, jl));
      vptdrydensInc_jlm1 = vptdrydensInc;
      vptdrydensInc = ((1.0 + constants::c_virtual * qView(jn, jl)) * ptIncView(jn, jl)
                       + ((1.0 + constants::c_virtual) * qIncView(jn, jl)
                       * ptView(jn, jl)
                       / (1.0 - qView(jn, jl))))
                      / (1.0 - qView(jn, jl));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensInc_intp_times_h_minus_hm1 = hl_minus_hm1 * vptdrydensInc +
                                              h_minus_hl * vptdrydensInc_jlm1;
      dryrhoIncView(jn, jl) = dryrhoView(jn, jl) * ((1.0 - constants::rd_over_cp)
                              * pIncView(jn, jl) / pView(jn, jl)
                              - vptdrydensInc_intp_times_h_minus_hm1
                              / vptdrydens_intp_times_h_minus_hm1);
    }
  }
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_new_dry_air_density_without_condensate_from_pressure_levels_ad(atlas::FieldSet
                                                                         & hatFlds,
                                                                         const atlas::FieldSet
                                                                         & stateFlds) {
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_ad()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto ptView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto qView = make_view<const double, 2>(stateFlds["specific_humidity"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  const auto dryrhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto ptHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto dryrhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;
  double vptdrydens;
  double vptdrydens_intp_times_h_minus_hm1;
  double vptdrydensHat;
  double vptdrydens_jlm1;
  double vptdrydensHat_jlm1;

  idx_t nlevels = hatFlds["dry_air_density_levels_minus_one"].levels();

  for (idx_t jn = 0; jn < dryrhoHatView.shape(0); ++jn) {
    for (idx_t jl = nlevels-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += dryrhoView(jn, jl) * (1.0 - constants::rd_over_cp) /
                          pView(jn, jl) * dryrhoHatView(jn, jl);
      vptdrydens = ptView(jn, jl) * (1.0 + constants::c_virtual * qView(jn, jl))
                   / (1.0 - qView(jn, jl));
      vptdrydens_jlm1 = ptView(jn, jl-1) * (1.0 + constants::c_virtual * qView(jn, jl-1))
                        / (1.0 - qView(jn, jl-1));
      vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                           hl_minus_hm1 * vptdrydens;
      vptdrydensHat = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * hl_minus_hm1
                      / vptdrydens_intp_times_h_minus_hm1;
      vptdrydensHat_jlm1 = - dryrhoView(jn, jl) * dryrhoHatView(jn, jl) * h_minus_hl
                           / vptdrydens_intp_times_h_minus_hm1;
      ptHatView(jn, jl) += (1.0 + constants::c_virtual * qView(jn, jl))
                           / (1.0 - qView(jn, jl))
                           * vptdrydensHat;
      qHatView(jn, jl) += (1.0 + constants::c_virtual) * ptView(jn, jl)
                          / ((1.0 - qView(jn, jl))
                          * (1.0 - qView(jn, jl)))
                          * vptdrydensHat;
      ptHatView(jn, jl-1) += (1.0 + constants::c_virtual * qView(jn, jl-1))
                             / (1.0 - qView(jn, jl-1))
                             * vptdrydensHat_jlm1;
      qHatView(jn, jl-1) += (1.0 + constants::c_virtual) * ptView(jn, jl-1)
                            / ((1.0 - qView(jn, jl-1))
                            * (1.0 - qView(jn, jl-1)))
                            * vptdrydensHat_jlm1;
      dryrhoHatView(jn, jl) = 0.0;
    }
    pHatView(jn, 0) += dryrhoView(jn, 0) * ((1.0 - constants::rd_over_cp) /
                       pView(jn, 0)) * dryrhoHatView(jn, 0);
    vptdrydens = ptView(jn, 0) * (1.0 + constants::c_virtual * qView(jn, 0))
                 / (1.0 - qView(jn, 0));
    vptdrydensHat = - dryrhoView(jn, 0) * dryrhoHatView(jn, 0)
                    / vptdrydens;
    ptHatView(jn, 0) += (1.0 + constants::c_virtual * qView(jn, 0))
                        / (1.0 - qView(jn, 0)) * vptdrydensHat;
    qHatView(jn, 0) += (1.0 + constants::c_virtual) * ptView(jn, 0)
                       / ((1.0 - qView(jn, 0))
                       * (1.0 - qView(jn, 0))) * vptdrydensHat;
    dryrhoHatView(jn, 0) = 0.0;
  }
  oops::Log::trace() <<
            "[eval_new_dry_air_density_without_condensate_from_pressure_levels_ad()] ... exit"
                     << std::endl;
}

}  // namespace mo
