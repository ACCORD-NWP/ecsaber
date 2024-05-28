/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/Increment.h"
#include "src/IncrModCtlVec.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// Covariance class

class Covariance : public util::Printable,
                   private eckit::NonCopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname()
    {return "quench::Covariance";}

/// OOPS interface

  Covariance(const Geometry &,
             const Variables &,
             const eckit::Configuration &,
             const State &)
    {}
  ~Covariance()
    {}

  void multiply(const Increment & dxi,
                Increment & dxo) const
    {dxo = dxi;}
  void inverseMultiply(const Increment & dxi,
                       Increment & dxo) const
    {dxo = dxi;}
  void randomize(Increment & dxo) const
    {dxo.random();}

 private:
  void print(std::ostream & os) const {os << "Covariance";}

/// ECSABER-specific definitions
#include "src/CovarianceECDef.h"
};

// -----------------------------------------------------------------------------

}  // namespace quench
