/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "vader/vader.h"
// Recipe headers
#include "recipes/AirDensityLevelsMinusOne.h"
#include "recipes/AirPotentialTemperature.h"
#include "recipes/AirPressureThickness.h"
#include "recipes/AirTemperature.h"
#include "recipes/AirVirtualTemperature.h"
#include "recipes/CloudIceMixingRatio.h"
#include "recipes/CloudLiquidMixingRatio.h"
#include "recipes/DryAirDensity.h"
#include "recipes/DryAirDensityLevelsMinusOne.h"
#include "recipes/HydrostaticExnerLevels.h"
#include "recipes/ParticulateMatter2p5.h"
#include "recipes/RainMixingRatio.h"
#include "recipes/TotalMixingRatio.h"
#include "recipes/TotalWater.h"
#include "recipes/uwind_at_10m.h"
#include "recipes/VirtualPotentialTemperature.h"
#include "recipes/vwind_at_10m.h"
#include "recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater.h"
#include "recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m.h"

namespace vader {

const cookbookConfigType Vader::defaultCookbookDefinition = {
        // Default VADER cookbook definition
        {oops::Variable{"air_temperature"},
                                     {AirTemperature_A::Name, AirTemperature_B::Name}},
        {oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"},
                                     {CloudIceMixingRatio_A::Name}},
        {oops::Variable{"cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"},
                                     {CloudLiquidMixingRatio_A::Name}},
        {oops::Variable{"dry_air_density"},
                                     {DryAirDensity_A::Name}},
        {oops::Variable{"dry_air_density_levels_minus_one"},
                                     {DryAirDensityLevelsMinusOne_A::Name}},
        {oops::Variable{"air_density_levels_minus_one"},
                                     {AirDensityLevelsMinusOne_A::Name,
                                      AirDensityLevelsMinusOne_B::Name}},
        {oops::Variable{"hydrostatic_exner_levels"},
                                     {HydrostaticExnerLevels_A::Name}},
        {oops::Variable{"potential_temperature"},
                                     {AirPotentialTemperature_A::Name}},
        {oops::Variable{"qrain"},    {RainMixingRatio_A::Name}},
        {oops::Variable{"m_t"},      {TotalMixingRatio_A::Name}},
        {oops::Variable{"qt"},       {TotalWater_A::Name}},
        {oops::Variable{"specific_humidity"},
                                     {WaterVaporMixingRatioWrtMoistAirAndCondensedWater_A::Name}},
        {oops::Variable{"specific_humidity_at_two_meters_above_surface"},
                                     {WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m_A::Name}},
        {oops::Variable{"uwind_at_10m"},
                                     {uwind_at_10m_A::Name}},
        {oops::Variable{"virtual_potential_temperature"},
                                     {VirtualPotentialTemperature_B::Name,
                                      VirtualPotentialTemperature_A::Name}},
        {oops::Variable{"virtual_temperature"},
                                     {AirVirtualTemperature_A::Name}},
        {oops::Variable{"vwind_at_10m"},
                                     {vwind_at_10m_A::Name}},
        {oops::Variable{"mass_density_of_particulate_matter_2p5_in_air"},
                                     {ParticulateMatter2p5_A::Name,
                                      ParticulateMatter2p5_B::Name}}
    };

}  // namespace vader
