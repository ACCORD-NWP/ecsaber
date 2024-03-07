/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field/for_each.h"

#include "mo/constants.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/functions.h"
#include "mo/lookups.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_sat_vapour_pressure_nl(const std::string & svp_file, atlas::FieldSet & fields)
{
  oops::Log::trace() << "[eval_sat_vapour_pressure_nl()] starting ..." << std::endl;

  // normalised T lambda function
  // which enforces upper and lower bounds on T
  auto normalisedT = [](const double tVal) {
    double t = (tVal - constants::TLoBound)/constants::Tinc;
    double t1 = (t < 0.0 ? 0.0 : t);
    double t2 = (t1 >= static_cast<double>(constants::svpLookUpLength - 1)) ?
                       static_cast<double>(constants::svpLookUpLength - 1) : t1;
    return t2;
  };

  // lambda function returning the index of a given temperature value;
  // avoids returning last index in table
  auto index = [](const double normalisedTVal) {
    std::size_t i = static_cast<std::size_t>(normalisedTVal);
    return (i == constants::svpLookUpLength - 1 ? constants::svpLookUpLength - 2 : i);
  };

  // weight lamba function
  auto weight = [index](const double normalisedTVal) {
    std::size_t i = index(normalisedTVal);
    return ( normalisedTVal -  static_cast<double>(i) );
  };

  // lambda function to interpolate values from table to given temperature
  auto interp = [](const double weight, const double v1, const double v2) {
     return (weight * v2 + (1 - weight) * v1);
  };

  // lambda function to determine the gradient of variable in table
  // with respect to temperature
  auto normalisedDiff = [](const double v1, const double v2) {
     return (v2 - v1)/ constants::Tinc;
  };

  double w;  // weight
  std::size_t indx;  // index

  const auto tView  = make_view<const double, 2>(fields["air_temperature"]);

  const std::vector<std::string> fnames {"svp", "dlsvpdT"};
  int ival = 0;  // set ival = 2 to get svp wrt water
  for (auto & ef : fnames) {
    if (fields.has(ef)) {
      // check this recipe to calculate svp is correct
      const auto& lookup = lookups::allSvpLookupTable[ival];
      ++ival;

      atlas::field::for_each_value(atlas::execution::par_unseq,
                                   fields["air_temperature"],
                                   fields[ef],
                                   [&](const double T,
                                       double& svp) {
        indx = index(normalisedT(T));
        w = weight(normalisedT(T));
        svp = interp(w, lookup.at(indx), lookup.at(indx+1));
      });
    }
  }

  oops::Log::trace() << "[eval_sat_vapour_pressure_nl()] ... exit" << std::endl;
}

}  // namespace mo

