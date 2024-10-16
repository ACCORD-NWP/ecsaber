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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/PressureParameters.h"

namespace oops {
  namespace patch {class Variables;}
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------
class HydrostaticExnerParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<GpToHpCovarianceParameters>
    hydrostaticexnerParams{"covariance data", this};
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables({
    "air_pressure_levels",
    "exner_levels_minus_one",
    "geostrophic_pressure_levels_minus_one",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels",
    "unbalanced_pressure_levels_minus_one"});}

  oops::JediVariables activeInnerVars(const oops::JediVariables& outerVars) const override {
    oops::JediVariables vars({"geostrophic_pressure_levels_minus_one",
                          "hydrostatic_exner_levels",
                          "hydrostatic_pressure_levels",
                          "unbalanced_pressure_levels_minus_one",
                         });
    const int modelLevels = outerVars.getLevels("exner_levels_minus_one");
    vars.addMetaData("geostrophic_pressure_levels_minus_one", "levels", modelLevels);
    vars.addMetaData("hydrostatic_exner_levels", "levels", modelLevels + 1);
    vars.addMetaData("hydrostatic_pressure_levels", "levels", modelLevels + 1);
    vars.addMetaData("unbalanced_pressure_levels_minus_one", "levels", modelLevels);
    return vars;
  }

  oops::JediVariables activeOuterVars(const oops::JediVariables& outerVars) const override {
    oops::JediVariables vars({"air_pressure_levels",
                          "hydrostatic_exner_levels",
                          "hydrostatic_pressure_levels",
                          "exner_levels_minus_one",
                         });
    vars.intersection(outerVars);

    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------
/// \brief  This saber block is here to do 3 jobs:
///         1) the vertical regression on geostrophic pressure
///         2) summing the result with unbalanced pressure to
///            create hydrostatic_pressure
///         3) converting hydrostatic pressure to exner pressure.
// TO DO: Marek - remove saber block when results identical to 3 new saber blocks
class HydrostaticExner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydrostaticExner";}

  typedef HydrostaticExnerParameters Parameters_;

  HydrostaticExner(const oops::GeometryData &,
                   const oops::JediVariables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~HydrostaticExner();

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
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
