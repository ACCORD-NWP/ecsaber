/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/generic/AtlasInterpolator.h"
#include "src/Covariance.h"
#include "src/Geometry.h"
#include "src/HorizScaleDecomposition.h"
#include "src/Increment.h"
#include "src/IncrEnsCtlVec.h"
#include "src/IncrModCtlVec.h"
#include "src/Interpolator.h"
#include "src/LinearVariableChange.h"
#include "src/LocalizationMatrix.h"
#include "src/Model.h"
#include "src/ModelAuxControl.h"
#include "src/ModelAuxControlEstimator.h"
#include "src/ModelAuxCovariance.h"
#include "src/ModelAuxCtlVec.h"
#include "src/ModelAuxIncrement.h"
#include "src/ModelData.h"
#include "src/State.h"
#include "src/VariableChange.h"
#include "src/Variables.h"


namespace oops {
class AtlasInterpolator;
}  // namespace oops

namespace quench {

struct Traits {
  static std::string name()
    {return "quench";}
  static std::string nameCovar()
    {return "quenchCovariance";}

  typedef quench::Covariance               Covariance;
  typedef quench::Geometry                 Geometry;
  typedef quench::HorizScaleDecomposition  HorizScaleDecomposition;
  typedef quench::Increment                Increment;
  typedef quench::IncrEnsCtlVec            IncrEnsCtlVec;
  typedef quench::IncrModCtlVec            IncrModCtlVec;
  typedef quench::Interpolator             Interpolator;
  typedef quench::LinearVariableChange     LinearVariableChange;
  typedef quench::LocalizationMatrix       LocalizationMatrix;
  typedef quench::Model                    Model;
  typedef quench::ModelAuxControl          ModelAuxControl;
  typedef quench::ModelAuxControlEstimator ModelAuxControlEstimator;
  typedef quench::ModelAuxCovariance       ModelAuxCovariance;
  typedef quench::ModelAuxCtlVec           ModelAuxCtlVec;
  typedef quench::ModelAuxIncrement        ModelAuxIncrement;
  typedef quench::ModelData                ModelData;
  typedef quench::State                    State;
  typedef quench::VariableChange           VariableChange;
  typedef quench::Variables                Variables;
};

}  // namespace quench
