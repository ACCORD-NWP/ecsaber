/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Increment;
  class IncrModCtlVec;
  class State;
  class Variables;

// -----------------------------------------------------------------------------
/// Covariance class

class Covariance : public util::Printable,
                   private boost::noncopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname()
    {return "quench::Covariance";}

  // Constructor/destructor
  Covariance(const Geometry &,
             const Variables &,
             const eckit::Configuration &,
             const State &)
    {}
  ~Covariance()
    {}

  // Linearize
  void linearize(const State &,
                 const Geometry &)
    {}

  // Multiply and inverse multiply
  void multiply(const Increment &,
                Increment &) const
    {throw eckit::NotImplemented(Here());}
  void inverseMultiply(const Increment &,
                       Increment &) const
    {throw eckit::NotImplemented(Here());}

  // Square-root and adjoint
  void multiplySqrt(const IncrModCtlVec &,
                    Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrtTrans(const Increment &,
                         IncrModCtlVec &) const
    {throw eckit::NotImplemented(Here());}

  // Randomization
  void randomize(Increment & dxo) const
    {throw eckit::NotImplemented(Here());}

 private:
  // Print
  void print(std::ostream & os) const
    {os << "Covariance";}
};

// -----------------------------------------------------------------------------

}  // namespace quench
