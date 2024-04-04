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
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

bool evalSatSpecificHumidity(atlas::FieldSet & stateFlds)
{
  oops::Log::trace() << "[evalSatSpecificHumidity()] starting ..." << std::endl;

  const auto & ghost = stateFlds["qsat"].functionspace().ghost();
  atlas::field::for_each_value_masked(atlas::execution::par_unseq,
                                      ghost,
                                      stateFlds["air_pressure"],
                                      stateFlds["svp"],
                                      stateFlds["air_temperature"],
                                      stateFlds["qsat"],
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

  stateFlds["qsat"].set_dirty();

  oops::Log::trace() << "[evalSatSpecificHumidity()] ... exit" << std::endl;

  return true;
}

}  // namespace mo
