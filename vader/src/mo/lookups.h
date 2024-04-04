/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mo/constants.h"
#include "mo/dlsvp_lookup.h"
#include "mo/dlsvpW_lookup.h"
#include "mo/svp_lookup.h"
#include "mo/svpW_lookup.h"

namespace mo {
namespace lookups {

static constexpr std::array<std::array<double, constants::svpLookUpLength>, 4> allSvpLookupTable {
  svpLookupTable, dlsvpLookupTable,
  svpWLookupTable, dlsvpWLookupTable
};

}  // namespace lookups
}  // namespace mo
