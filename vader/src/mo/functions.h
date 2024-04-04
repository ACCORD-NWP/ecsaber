/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace mo {
namespace functions {
//--
// ++ Atlas Function Spaces ++

/// \brief procedure to call a functor for a given concrete implementation
///        of a function space type
template<typename Functor>
void executeFunc(const atlas::FunctionSpace & fspace, const Functor & functor) {
  if (atlas::functionspace::NodeColumns(fspace)) {
    functor(atlas::functionspace::CubedSphereNodeColumns(fspace));
  } else if (atlas::functionspace::CellColumns(fspace)) {
    functor(atlas::functionspace::CubedSphereCellColumns(fspace));
  } else if (atlas::functionspace::StructuredColumns(fspace)) {
    functor(atlas::functionspace::StructuredColumns(fspace));
  } else {
    oops::Log::error() << "ERROR - a functor call failed "
                          "(function space type not allowed)" << std::endl;
    throw std::runtime_error("a functor call failed");
  }
}

/// \brief wrapper for 'parallel_for'
template<typename Functor>
void parallelFor(const atlas::FunctionSpace & fspace,
                 const Functor& functor,
                 const atlas::util::Config& conf = atlas::util::Config()) {
  executeFunc(fspace, [&](const auto& fspace){fspace.parallel_for(conf, functor);});
}


//--
// ++ I/O processing ++

/// \details This extracts the scaling coefficients that are applied to Cleff and Cfeff
///          to generate the qcl and qcf increments in the moisture incrementing operator (MIO)
///          The string s can be "qcl_coef" or "qcf_coef"
Eigen::MatrixXd createMIOCoeff(const std::string mioFileName,
                               const std::string s);

extern "C" {
  void umGetLookUp2D_f90(
    const int &,
    const char *,
    const int &,
    const char *,
    const int &,
    const int &,
    double &);
}  // extern "C"

}  // namespace functions
}  // namespace mo
