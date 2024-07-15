/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/assimilation/instantiateMinFactory.h"

#include "saber/oops/DRPLanczosEVILMinimizer.h"
#include "saber/oops/PLanczosEVILMinimizer.h"
#include "saber/oops/RPLanczosEVILMinimizer.h"
#include "saber/oops/SQRTBPLanczosEVILMinimizer.h"
#include "saber/oops/SQRTPLanczosEVILMinimizer.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateMinFactory() {
  static oops::MinMaker<MODEL, DRPLanczosEVILMinimizer<MODEL> >
    makerDRPLanczosEVIL_("DRPLanczosEVIL");
  static oops::MinMaker<MODEL, PLanczosEVILMinimizer<MODEL> > makerPLanczosEVIL_("PLanczosEVIL");
  static oops::MinMaker<MODEL, RPLanczosEVILMinimizer<MODEL> > makerRPLanczosEVIL_("RPLanczosEVIL");
  static oops::MinMaker<MODEL, SQRTBPLanczosEVILMinimizer<MODEL> >
    makerSQRTBPLanczosEVIL_("SQRTBPLanczosEVIL");
  static oops::MinMaker<MODEL, SQRTPLanczosEVILMinimizer<MODEL> >
    makerSQRTPLanczosEVIL_("SQRTPLanczosEVIL");
}

// -----------------------------------------------------------------------------

}  // namespace saber
