/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace quench {

class Covariance;
class Geometry;
class GeoVaLs;
class HorizScaleDecomposition;
class Increment;
class IncrEnsCtlVec;
class IncrModCtlVec;
class Interpolator;
class LocalizationMatrix;
class Locations;
class Model;
class ModelAuxControl;
class ModelAuxControlEstimator;
class ModelAuxCovariance;
class ModelAuxCtlVec;
class ModelAuxIncrement;
class ObsSpace;
class ObsVector;
class State;
class Variables;

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
  using ObsVector = quench::ObsVector;
  using State = quench::State;
  using Variables = quench::Variables;
};

}  // namespace quench
