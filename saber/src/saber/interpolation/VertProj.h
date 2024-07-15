/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

// Note that this is a saber block to demonstrate proof of concept
// When we have a fully working vertical mode to model level saber block
// we can get rid of this.
namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------
class VertProjParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VertProjParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<oops::JediVariables> activeVariables{"active variables", this};
  oops::RequiredParameter<atlas::idx_t> innerVerticalLevels{"inner vertical levels",
    "inner number of vertical levels", this};
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables();}
};

// -----------------------------------------------------------------------------

class VertProj : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::VertProj";}

  typedef VertProjParameters Parameters_;

  VertProj(const oops::GeometryData &,
           const oops::JediVariables &,
           const eckit::Configuration &,
           const Parameters_ &,
           const oops::FieldSet3D &,
           const oops::FieldSet3D &);
  virtual ~VertProj() = default;

  const oops::GeometryData & innerGeometryData()
    const override {return outerGeometryData_;}
  const oops::JediVariables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  const oops::GeometryData & outerGeometryData_;
  oops::JediVariables outerVars_;
  oops::JediVariables activeVars_;
  oops::JediVariables innerVars_;
};

}  // namespace interpolation
}  // namespace saber
