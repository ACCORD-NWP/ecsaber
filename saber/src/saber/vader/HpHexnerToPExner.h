/*
 * (C) Crown Copyright 2023 Met Office
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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace oops {
  class JediVariables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

class HpHexnerToPExnerParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HpHexnerToPExnerParameters, SaberBlockParametersBase)

 public:
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables({
    std::vector<std::string>{
    "air_pressure_levels",
    "exner_levels_minus_one",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels"}});}

  oops::JediVariables activeInnerVars(const oops::JediVariables& outerVars) const override {
    const int modelLevels = outerVars["air_pressure_levels"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::JediVariables vars;
    vars.push_back({"hydrostatic_exner_levels", conf});
    vars.push_back({"hydrostatic_pressure_levels", conf});
    return vars;
  }

  oops::JediVariables activeOuterVars(const oops::JediVariables& outerVars) const override {
    oops::JediVariables vars({outerVars["air_pressure_levels"],
                          outerVars["exner_levels_minus_one"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------
/// \brief This saber block is here to copy hydrostatic_exner_levels into
///        exner_levels_minus_one and hydrostatic_pressure_levels into
///        air_pressure_levels

class HpHexnerToPExner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HpHexnerToPExner";}

  typedef HpHexnerToPExnerParameters Parameters_;

  HpHexnerToPExner(const oops::GeometryData &,
                   const oops::JediVariables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~HpHexnerToPExner();

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
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
