/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

class AirTemperature_A_Parameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperature_A_Parameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

class AirTemperature_B_Parameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperature_B_Parameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

class AirTemperature_C_Parameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperature_C_Parameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

// ------------------------------------------------------------------------------------------------
/* CCPP Names for (hopefully near) future:
 *  \brief AirTemperature_A class defines a recipe for air_temperature_at_interface
 *         from air_potential_temperature_at_interface and
 *         dimensionless_exner_function_at_interface
 *
 *         NL air_temperature_at_interface = air_potential_temperature_at_interface
 *                                              * dimensionless_exner_function_at_interface
 */
/*! \brief AirTemperature_A class defines a recipe for air_temperature
 *         from potential_temperature and exner
 *
 *         NL air_temperature = potential_temperature * exner
 */
class AirTemperature_A : public RecipeBase {
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef AirTemperature_A_Parameters Parameters_;

    AirTemperature_A(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};

// ------------------------------------------------------------------------------------------------
/*! \brief AirTemperature_B class defines a recipe for temperature from virtual temperature
 *         and specific humidity.
 */
class AirTemperature_B : public RecipeBase {
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef AirTemperature_B_Parameters Parameters_;

    AirTemperature_B(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

// ------------------------------------------------------------------------------------------------
/*! \brief AirTemperature_C class defines a recipe for temperature from air pressure and air potential
 *         temperature.
 *
 *  \detail This recipe is a special case for WRF. It calculates air potential temperature using 
 *          perturbation air potential temperature (T in wrfout) and base air potential temperature 
 *          (pt_base, default 290K) specified by a yaml key. Then multiply exner function derived 
 *          by air pressure to get air temperature. It does not provide TL/AD algorithms.
 */
class AirTemperature_C : public RecipeBase {
 public:
    static const char Name[];
    static const oops::JediVariables Ingredients;

    typedef AirTemperature_C_Parameters Parameters_;

    AirTemperature_C(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    oops::Variable product() const override;
    oops::JediVariables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

}  // namespace vader

