/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>

#include "quench/Covariance.h"
#include "quench/Geometry.h"
#include "quench/GeoVaLs.h"
#include "quench/HorizScaleDecomposition.h"
#include "quench/Increment.h"
#include "quench/IncrEnsCtlVec.h"
#include "quench/IncrModCtlVec.h"
#include "quench/Interpolator.h"
#include "quench/LocalizationMatrix.h"
#include "quench/Locations.h"
#include "quench/Model.h"
#include "quench/ModelAuxControl.h"
#include "quench/ModelAuxControlEstimator.h"
#include "quench/ModelAuxCovariance.h"
#include "quench/ModelAuxCtlVec.h"
#include "quench/ModelAuxIncrement.h"
#include "quench/ObsSpace.h"
#include "quench/ObsVec.h"
#include "quench/State.h"
#include "quench/Variables.h"

namespace quench {

// -----------------------------------------------------------------------------

struct Traits {
  static std::string name()
    {return "quench";}
  static std::string nameCovar()
    {return "quenchCovariance";}

  using Covariance = quench::Covariance;
  using Geometry = quench::Geometry;
  using GeoVaLs = quench::GeoVaLs;
  using HorizScaleDecomposition = quench::HorizScaleDecomposition;
  using Increment = quench::Increment;
  using IncrEnsCtlVec = quench::IncrEnsCtlVec;
  using IncrModCtlVec = quench::IncrModCtlVec;
  using Interpolator = quench::Interpolator;
  using LocalizationMatrix = quench::LocalizationMatrix;
  using Locations = quench::Locations;
  using Model = quench::Model;
  using ModelAuxControl = quench::ModelAuxControl;
  using ModelAuxControlEstimator = quench::ModelAuxControlEstimator;
  using ModelAuxCovariance = quench::ModelAuxCovariance;
  using ModelAuxCtlVec = quench::ModelAuxCtlVec;
  using ModelAuxIncrement = quench::ModelAuxIncrement;
  using ObsSpace = quench::ObsSpace;
  using ObsVector = quench::ObsVec;
  using State = quench::State;
  using Variables = quench::Variables;
};

// -----------------------------------------------------------------------------

}  // namespace quench
