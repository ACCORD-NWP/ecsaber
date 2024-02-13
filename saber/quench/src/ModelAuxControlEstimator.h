/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Model;
  class ModelAuxControl;
  class State;

// -----------------------------------------------------------------------------
/// ModelAuxControlEstimator class

class ModelAuxControlEstimator : public util::Printable,
                             private eckit::NonCopyable,
                             private util::ObjectCounter<ModelAuxControlEstimator> {
 public:
  static const std::string classname()
    {return "quench::ModelAuxControlEstimator";}

/// Constructors, destructors
  ModelAuxControlEstimator(const Geometry &,
                           const Model &,
                           const eckit::Configuration &)
    {}
  ~ModelAuxControlEstimator()
    {}

/// Methods
  void estimate(const State &,
                ModelAuxControl &) const
    {}
  void estimate(const State &,
                const State &,
                ModelAuxControl &) const
    {}

 private:
  void print(std::ostream & os) const
    {}
};

// -----------------------------------------------------------------------------

}  // namespace quench
