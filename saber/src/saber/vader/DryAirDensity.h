/*
 * (C) Crown Copyright 2022-2024 Met Office
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

class DryAirDensityParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensityParameters, SaberBlockParametersBase)

 public:
  oops::JediVariables mandatoryActiveVars() const override {return oops::JediVariables({
    std::vector<std::string>{
    "dry_air_density_levels_minus_one",
    "air_pressure_levels",
    "potential_temperature",
    "specific_humidity",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "mass_content_of_cloud_ice_in_atmosphere_layer"}});}

  const oops::JediVariables mandatoryStateVars() const override {return oops::JediVariables({
    "dry_air_density_levels_minus_one",
    "air_pressure_levels_minus_one",
    "potential_temperature",
    "specific_humidity",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "mass_content_of_cloud_ice_in_atmosphere_layer"});}

  oops::JediVariables activeInnerVars(const oops::JediVariables& outerVars) const override {
    const int modelLevels = outerVars["dry_air_density_levels_minus_one"].getLevels();
    oops::JediVariables vars;
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels + 1);
    vars.push_back({"air_pressure_levels", conf});
    conf.set("levels", modelLevels);
    vars.push_back({"potential_temperature", conf});
    vars.push_back({"specific_humidity", conf});
    vars.push_back({"mass_content_of_cloud_liquid_water_in_atmosphere_layer", conf});
    vars.push_back({"mass_content_of_cloud_ice_in_atmosphere_layer", conf});
    return vars;
  }

  oops::JediVariables activeOuterVars(const oops::JediVariables& outerVars) const override {
    oops::JediVariables vars({outerVars["dry_air_density_levels_minus_one"]});
    return vars;
  }

  oops::JediVariables intermediateTempVars(const oops::JediVariables& outerVars) const {
    if (outerVars.has("air_pressure_levels_minus_one")) {
      throw eckit::UserError("air_pressure_levels_minus_one is a "
                             "temporary variable of mo_dry_air_density "
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels = outerVars["dry_air_density_levels_minus_one"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::JediVariables tempVars;
    tempVars.push_back({"air_pressure_levels_minus_one", conf});
    return tempVars;
  }
};

// -----------------------------------------------------------------------------

class DryAirDensity : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::DryAirDensity";}

  typedef DryAirDensityParameters Parameters_;

  DryAirDensity(const oops::GeometryData &,
                const oops::JediVariables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const oops::FieldSet3D &,
                const oops::FieldSet3D &);
  virtual ~DryAirDensity();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::JediVariables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::JediVariables innerVars_;
  const oops::JediVariables activeOuterVars_;
  const oops::JediVariables innerOnlyVars_;
  const oops::JediVariables intermediateTempVars_;
  const oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
