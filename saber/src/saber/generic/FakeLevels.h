/*
 * (C) Copyright 2024 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------
class FakeLevelsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(FakeLevelsParameters, SaberBlockParametersBase)

 public:
  // Fake levels
  oops::RequiredParameter<std::vector<double>> fakeLevels{"fake levels", this};

  // Input model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputModelFilesConf{
    "input model files", this};

  // Scalar vertical support
  oops::OptionalParameter<double> rvFromYaml{"vertical length-scale", this};

  // Mandatory variables
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class FakeLevels : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::generic::FakeLevels";}

  typedef FakeLevelsParameters Parameters_;

  FakeLevels(const oops::GeometryData &,
             const oops::patch::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const oops::FieldSet3D &,
             const oops::FieldSet3D &);
  virtual ~FakeLevels() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return gdata_;}
  const oops::patch::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;

 private:
  oops::patch::Variables createInnerVars(const atlas::idx_t &,
                                  const oops::patch::Variables &,
                                  const oops::patch::Variables &) const;
  void print(std::ostream &) const override;

  Parameters_ params_;
  size_t nz_;
  const util::DateTime validTime_;
  const oops::GeometryData & gdata_;
  const eckit::mpi::Comm & comm_;
  oops::patch::Variables outerVars_;
  oops::patch::Variables activeVars_;
  const std::string suffix_;
  oops::patch::Variables innerVars_;
  std::vector<double> fakeLevels_;
  std::unique_ptr<oops::FieldSet3D> weights_;
  eckit::LocalConfiguration fieldsMetaData_;
};

}  // namespace generic
}  // namespace saber
