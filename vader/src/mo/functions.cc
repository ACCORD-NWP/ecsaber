/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Core>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/constants.h"
#include "mo/functions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {
namespace functions {

Eigen::MatrixXd createMIOCoeff(const std::string mioFileName,
                               const std::string s) {
  oops::Log::trace() << "[createMIOCoeff] starting ..." << std::endl;
  Eigen::MatrixXd mioCoeff(static_cast<std::size_t>(constants::mioLevs),
                           static_cast<std::size_t>(constants::mioBins));

  std::vector<double> valuesvec(constants::mioLookUpLength, 0);

  umGetLookUp2D_f90(static_cast<int>(mioFileName.size()),
                    mioFileName.c_str(),
                    static_cast<int>(s.size()),
                    s.c_str(),
                    static_cast<int>(constants::mioBins),
                    static_cast<int>(constants::mioLevs),
                    valuesvec[0]);

  for (int j = 0; j < static_cast<int>(constants::mioLevs); ++j) {
    for (int i = 0; i < static_cast<int>(constants::mioBins); ++i) {
      // Fortran returns column major order, but C++ needs row major
      mioCoeff(j, i) = valuesvec[i * constants::mioLevs+j];
    }
  }
  oops::Log::trace() << "[createMIOCoeff] ... exit" << std::endl;
  return mioCoeff;
}

}  // namespace functions
}  // namespace mo
