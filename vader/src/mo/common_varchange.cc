/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/for_each.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "mo/common_varchange.h"
#include "mo/constants.h"
#include "mo/functions.h"
#include "mo/lookups.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

bool evalSatSpecificHumidity(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalSatSpecificHumidity()] starting ..." << std::endl;

  atlas::field::for_each_value(atlas::execution::par_unseq,
                               fields["air_pressure"],
                               fields["svp"],
                               fields["air_temperature"],
                               fields["qsat"],
                               [&](const double pBar,
                                   const double svp,
                                   const double T,
                                   double& qsat) {
    // Compute the factor that converts from sat vapour pressure in a
    // pure water system to sat vapour pressure in air, FSUBW. This formula
    // is taken from equation A4.7 of Adrian Gill's book: Atmosphere-Ocean
    // Dynamics. Note that his formula works in terms of pressure in mb and
    // temperature in Celsius, so conversion of units leads to the slightly
    // different equation used here.
    const double fsubw = 1.0 + 1.0E-8 * pBar * (4.5 + 6.0e-4
                          * (T - constants::zerodegc)
                          * (T - constants::zerodegc));

    // Multiply by FSUBW to convert to saturated vapour pressure in air
    // (equation A4.6 of Adrian Gill's book)
    const double svp_times_fsubw = fsubw * svp;

    // Now form the accurate expression for QS, which is a rearranged
    // version of equation A4.3 of Gill's book.
    // Note that at very low pressures we apply a fix, to prevent a
    // singularity (Qsat tends to 1.0 kg/kg).
    qsat = constants::rd_over_rv * svp_times_fsubw /
          (std::max(pBar, svp_times_fsubw) -
          (1.0 - constants::rd_over_rv) * svp_times_fsubw);
  });

  oops::Log::trace() << "[evalSatSpecificHumidity()] ... exit" << std::endl;

  return true;
}

bool evalAirPressureLevels(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalAirPressureLevels()] starting ..." << std::endl;

  const auto ds_elmo = make_view<const double, 2>(fields["exner_levels_minus_one"]);
  const auto ds_plmo = make_view<const double, 2>(fields["air_pressure_levels_minus_one"]);
  const auto ds_t = make_view<const double, 2>(fields["potential_temperature"]);
  const auto ds_hl = make_view<const double, 2>(fields["height_levels"]);
  auto ds_pl = make_view<double, 2>(fields["air_pressure_levels"]);

  idx_t levels(fields["air_pressure_levels"].shape(1));
  for (idx_t jn = 0; jn < fields["air_pressure_levels"].shape(0); ++jn) {
    for (idx_t jl = 0; jl < levels - 1; ++jl) {
      ds_pl(jn, jl) = ds_plmo(jn, jl);
    }

    // Note that I am calculating the exner pressure above the top first and then
    // converting it to pressure
    // Note that strictly speaking we should be using virtual potential temperature here
    // but given that there should be no moisture at the top of the model, we should be
    // able to ignore that contribution.
    //
    // exner^k+1 = exner^k - g (height^k+1 - height^k) / (cp theta_v)
    //
    // pressure^k+1 = reference_pressure * (exner^k+1)**((1.0 / constants::rd_over_cp)
    //
    // where k is the model level index on half levels just below model top.

    ds_pl(jn, levels-1) =  constants::p_zero * pow(
      ds_elmo(jn, levels-2) - (constants::grav * (ds_hl(jn, levels-1) - ds_hl(jn, levels-2))) /
      (constants::cp * ds_t(jn, levels-2)), (1.0 / constants::rd_over_cp));

    ds_pl(jn, levels-1) = ds_pl(jn, levels-1) > 0.0 ? ds_pl(jn, levels-1) : constants::deps;
  }

  oops::Log::trace() << "[evalAirPressureLevels()] ... exit" << std::endl;

  return true;
}

}  // namespace mo
