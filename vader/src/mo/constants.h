/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <string>

namespace mo {
namespace constants {
  // Lower temperature bound in lookup table
  static constexpr double TLoBound = 183.15;

  // Upper temperature bound in lookup table
  static constexpr double THiBound = 338.15;

  // Temperature increment in lookup table
  static constexpr double Tinc = 0.1;

  // Round off coefficient for static cast
  static constexpr double deps = 1e-12;

  // Size of Look up table (1551)
  static constexpr int svpLookUpLength =
    static_cast<int>((THiBound - TLoBound + deps) / Tinc) + 1;

  static constexpr std::size_t mioBins = 21;
  static constexpr std::size_t mioLevs = 40;
  static constexpr std::size_t mioLookUpLength = mioLevs * mioBins;
  static constexpr double rHTBin = 5.0;     // done as a percentage!
  static constexpr double MinRhRef = 0.0;
  static constexpr double MaxRhRef = 150.0;
  static constexpr double rHTLastBinLowerLimit = (mioBins - 1) * rHTBin;  // in percentage
  static constexpr double TolMonotonicity = 1.0;
  static constexpr double effectiveRNegative = 5.0;

  // tolerance for avoiding division by zero in eval_mio_fields_nl
  static constexpr double tol = 1.0e-5;

  static constexpr double p_zero = 1.0e5;
  static constexpr double zerodegc = 273.15;             // conversion between degrees Celsius
                                                         // and Kelvin
  static constexpr double deg2rad        = M_PI / 180.;
  static constexpr double rad2deg        = 180. * M_1_PI;
  static constexpr double grav           = 9.80665e+0;
  static constexpr double t0c            = 2.7315e+2;    // temperature at zero celsius (K)
  static constexpr double ttp            = 2.7316e+2;    // temperature at h2o triple point (K)
  static constexpr double rd             = 2.8705e2;
  static constexpr double rv             = 4.6150e2;
  static constexpr double cp             = 1.0046e3;     // heat capacity at constant pressure
                                                         //      for air
  static constexpr double cv             = 7.1760e2;     // heat capacity at constant volume
                                                         //      for air
  static constexpr double rspec          = cp - cv;      // specific gas constant for dry air
  static constexpr double rspec_over_cp  = rspec / cp;
  static constexpr double pref           = 1.0e5;        // Reference pressure for calculating
                                                         //      exner
  static constexpr double rd_over_rv     = rd / rv;
  static constexpr double rd_over_cp     = rd / cp;
  static constexpr double cv_over_cp     = cv / cp;
  static constexpr double rv_over_rd     = rv / rd;
  static constexpr double rd_over_g      = rd / grav;
  static constexpr double g_over_rd      = grav / rd;
  static constexpr double mean_earth_rad = 6371.0;
  static constexpr double zero           = 0.0;
  static constexpr double quarter        = 0.25;
  static constexpr double half           = 0.5;
  static constexpr double one            = 1.0;
  static constexpr double two            = 2.0;
  static constexpr double four           = 4.0;
  static constexpr double five           = 5.0;
  static constexpr double ten            = 10.0;
  static constexpr double k_t            = 0.65;         // Thermal conductivity of water
  static constexpr double L_e            = 2.26e+06;     // Latent heat of vaporization
  static constexpr double eps            = 0.1;          // Albedo of sea water
  static constexpr double sig            = 5.67e-6;      // Stefan-Boltzmann constant
  static constexpr double alpha          = 2.7e-4;       // Water thermal expansion coefficient
  static constexpr double cw             = 0.015;        // Water specific heat
  static constexpr double v_w            = 0.8e-6;       // Water kinematic viscosity
  static constexpr double S_B            = 0.026;
  static constexpr double gr             = 9.81;
  static constexpr double Rou            = 1000.0;
  static constexpr double DU             = 21.4e-6;      // Dobson unit, kg O3/m**2
  static constexpr double Lclr           = 0.0065;       // constant lapse rate
  static constexpr double t2tv           = 0.608;        // constant lapse rate
  static constexpr double von_karman     = 0.41;         // Von Karman Constant
  static constexpr double es_w_0         = 611.2;        // saturation vapor pressure of water
                                                         //   at 0degC
  static constexpr double euzc_0         = 34.0;         // constant for estimating euphotic layer
  static constexpr double euzc_1         = -0.39;        // constant for estimating euphotic layer
  static constexpr double epsilon        = 0.62198;      // ratio of molecular weight of
                                                         //   water and dry air
  static constexpr double c_virtual      = (1.0 / epsilon) - 1.0;

  // International Civil Aviation Organization (ICAO) atmosphere.
  // https://en.wikipedia.org/wiki/International_Standard_Atmosphere#ICAO_Standard_Atmosphere
  // with gpm = height in geopotential metres
  static constexpr double icao_lapse_rate_l = 6.5E-03;   // Lapse rate for levels up
                                                         //      to 11,000 gpm
                                                         //      (lower boundary of the
                                                         //      isothermal layer) [K/gpm]
  static constexpr double icao_lapse_rate_u = -1.0E-03;  // Lapse rate for levels above
                                                         //      20,000 gpm
                                                         //      (upper boundary of the
                                                         //      isothermal layer) [K/gpm]
  static constexpr double icao_height_l     = 11000.0;   // Height of bottom of isothermal
                                                         //      layer [gpm]
  static constexpr double icao_height_u     = 20000.0;   // Height of top of isothermal
                                                         //      layer [gpm]
  static constexpr double icao_temp_surface = 288.15;    // Surface temperature [K]
  static constexpr double icao_temp_isothermal_layer = 216.65;  // Temperature of isothermal
                                                                //      layer [K]
  static constexpr double icao_pressure_surface = 1013.25;  // Assumed surface pressure [hPa]
  static constexpr double icao_pressure_l   = 226.32;    // Assumed pressure at 11,000 gpm [hPa]
  static constexpr double icao_pressure_u   = 54.7487;   // Assumed pressure at 20,000 gpm [hPa]

}  // namespace constants
}  // namespace mo
