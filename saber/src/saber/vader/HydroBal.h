/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace oops {
  class JediVariables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HydroBalParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydroBalParameters, SaberBlockParametersBase)

 public:
  oops::JediVariables mandatoryActiveVars() const override {
    return oops::JediVariables({
        std::vector<std::string>{
        "air_pressure_levels",
        "hydrostatic_exner_levels",
        "virtual_potential_temperature"}});
  }

  const oops::JediVariables mandatoryStateVars() const override {
    return oops::JediVariables({
        "air_pressure_levels",
        "hydrostatic_exner_levels",
        "virtual_potential_temperature",
        "height_levels"});
  }

  oops::JediVariables activeInnerVars(const oops::JediVariables& outerVars) const override {
    const int modelLevels = outerVars["virtual_potential_temperature"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels + 1);
    oops::JediVariables vars;
    vars.push_back({"air_pressure_levels", conf});
    vars.push_back({"hydrostatic_exner_levels", conf});
    return vars;
  }

  oops::JediVariables activeOuterVars(const oops::JediVariables& outerVars) const override {
    oops::JediVariables vars({outerVars["virtual_potential_temperature"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class HydroBal : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydroBal";}

  typedef HydroBalParameters Parameters_;

  HydroBal(const oops::GeometryData &,
           const oops::JediVariables &,
           const eckit::Configuration &,
           const Parameters_ &,
           const oops::FieldSet3D &,
           const oops::FieldSet3D &);
  virtual ~HydroBal();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::JediVariables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::JediVariables innerVars_;
  const oops::JediVariables activeOuterVars_;
  const oops::JediVariables innerOnlyVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
