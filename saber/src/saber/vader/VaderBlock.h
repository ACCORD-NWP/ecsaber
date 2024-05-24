/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

#include "vader/vader.h"

namespace saber {

// -----------------------------------------------------------------------------

class VaderBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VaderBlockParameters, SaberBlockParametersBase)

 public:
  oops::Parameter<vader::VaderParameters> vader{"vader", {}, this};
  oops::RequiredParameter<oops::JediVariables> innerVars{"inner variables", this};
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables();}
};

// -----------------------------------------------------------------------------

class VaderBlock : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::VaderBlock";}

  typedef VaderBlockParameters Parameters_;

  VaderBlock(const oops::GeometryData &,
             const oops::JediVariables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const oops::FieldSet3D &,
             const oops::FieldSet3D &);

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::JediVariables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  const oops::JediVariables outerVars_;
  const oops::GeometryData & innerGeometryData_;
  const oops::JediVariables innerVars_;
  vader::Vader vader_;
};

}  // namespace saber
