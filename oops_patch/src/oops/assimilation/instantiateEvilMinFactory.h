/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/assimilation/instantiateMinFactory.h"

#include "oops/assimilation/DRPLanczosEVILMinimizer.h"
#include "oops/assimilation/PLanczosEVILMinimizer.h"
#include "oops/assimilation/RPLanczosEVILMinimizer.h"
#include "oops/assimilation/SQRTBPLanczosEVILMinimizer.h"
#include "oops/assimilation/SQRTPLanczosEVILMinimizer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateEvilMinFactory() {
  static MinMaker<MODEL, DRPLanczosEVILMinimizer<MODEL> > makerDRPLanczosEVIL_("DRPLanczosEVIL");
  static MinMaker<MODEL, PLanczosEVILMinimizer<MODEL> > makerPLanczosEVIL_("PLanczosEVIL");
  static MinMaker<MODEL, RPLanczosEVILMinimizer<MODEL> > makerRPLanczosEVIL_("RPLanczosEVIL");
  static MinMaker<MODEL, SQRTBPLanczosEVILMinimizer<MODEL> >
    makerSQRTBPLanczosEVIL_("SQRTBPLanczosEVIL");
  static MinMaker<MODEL, SQRTPLanczosEVILMinimizer<MODEL> >
    makerSQRTPLanczosEVIL_("SQRTPLanczosEVIL");
}

// -----------------------------------------------------------------------------

}  // namespace oops
