/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>

#include "src/Covariance.h"
#include "src/Geometry.h"
#include "src/GeoVals.h"
#include "HorizScaleDecomposition.h"
#include "src/Increment.h"
#include "src/IncrEnsCtlVec.h"
#include "src/IncrModCtlVec.h"
#include "src/Interpolator.h"
#include "src/LocalizationMatrix.h"
#include "src/Locations.h"
#include "src/Model.h"
#include "src/ModelAuxControl.h"
#include "src/ModelAuxControlEstimator.h"
#include "src/ModelAuxCovariance.h"
#include "src/ModelAuxCtlVec.h"
#include "src/ModelAuxIncrement.h"
#include "src/ObsSpace.h"
#include "src/ObsVec.h"
#include "src/State.h"
#include "src/Variables.h"

namespace quench {

// -----------------------------------------------------------------------------

struct Traits {
  static std::string name()
    {return "quench";}
  static std::string nameCovar()
    {return "quenchCovariance";}

  using Geometry = quench::Geometry;

  using State = quench::State;
  using Model = quench::Model;
  using Increment = quench::Increment;
  using IncrEnsCtlVec = quench::IncrEnsCtlVec;
  using IncrModCtlVec = quench::IncrModCtlVec;
  using Covariance = quench::Covariance;
  using Variables = quench::Variables;

  using ModelAuxControl = quench::ModelAuxControl;
  using ModelAuxControlEstimator = quench::ModelAuxControlEstimator;
  using ModelAuxCovariance = quench::ModelAuxCovariance;
  using ModelAuxCtlVec = quench::ModelAuxCtlVec;
  using ModelAuxIncrement = quench::ModelAuxIncrement;

  using ObsSpace = quench::ObsSpace;
  using ObsVector = quench::ObsVec;

  using GeoVaLs = quench::GeoVaLs;
  using Locations = quench::Locations;

  using LocalizationMatrix = quench::LocalizationMatrix;
  using Interpolator = quench::Interpolator;
  using HorizScaleDecomposition = quench::HorizScaleDecomposition;
};

// -----------------------------------------------------------------------------

}  // namespace quench
