/*
 * (C) Crown Copyright 2022-2023 Met Office
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbCalibrationVertCovParameters>
    calibrationParams{"calibration", this};
  oops::OptionalParameter<spectralbReadParameters> readParams{"read", this};
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables();}
};

// -----------------------------------------------------------------------------

class SpectralCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralCovariance";}

  typedef SpectralCovarianceParameters Parameters_;

  SpectralCovariance(const oops::GeometryData &,
                     const oops::JediVariables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const oops::FieldSet3D &,
                     const oops::FieldSet3D &);

  virtual ~SpectralCovariance() = default;

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void write() const override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::JediVariables activeVars_;
  /// Vertical Spectral Covariances
  atlas::FieldSet spectralVerticalCovariances_;
  /// Geometry data
  const oops::GeometryData & geometryData_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
};

}  // namespace spectralb
}  // namespace saber
