/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/assimilation/instantiateMinFactory.h"

#include "oops/assimilation/DRPLanczosEVILMinimizer_cy46.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateEvilMinFactory() {
  static MinMaker<MODEL, DRPLanczosEVILMinimizer<MODEL> > makerDRPLanczosEVIL_("DRPLanczosEVIL");
}

// -----------------------------------------------------------------------------

}  // namespace oops
