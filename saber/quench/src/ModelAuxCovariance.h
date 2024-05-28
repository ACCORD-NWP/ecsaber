/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class ModelAuxControl;
  class ModelAuxCtlVec;
  class ModelAuxIncrement;
  class Geometry;

// -----------------------------------------------------------------------------
/// ModelAuxCovariance class

class ModelAuxCovariance : public util::Printable,
                           private eckit::NonCopyable,
                           private util::ObjectCounter<ModelAuxCovariance> {
 public:
  static const std::string classname()
    {return "quench::ModelAuxCovariance";}

/// OOPS interface

// Constructor/destructor
  ModelAuxCovariance(const eckit::Configuration &,
                     const Geometry &)
    {}
  ~ModelAuxCovariance()
    {}

// Linear algebra operators
  void linearize(const ModelAuxControl &,
                 const Geometry &)
    {}
  void multiply(const ModelAuxIncrement &,
                ModelAuxIncrement &) const
    {}
  void inverseMultiply(const ModelAuxIncrement &,
                       ModelAuxIncrement &) const
    {}
  void multiplySqrt(const ModelAuxCtlVec &,
                    ModelAuxIncrement &) const
    {}
  void multiplySqrtTrans(const ModelAuxIncrement &,
                         ModelAuxCtlVec &) const
    {}
  void randomize(ModelAuxIncrement &) const
    {}

// Return configuration
  const eckit::Configuration & config() const
    {return conf_;}

 private:
  void print(std::ostream & os) const
    {}

  const eckit::LocalConfiguration conf_ = eckit::LocalConfiguration();
};

// -----------------------------------------------------------------------------

}  // namespace quench
