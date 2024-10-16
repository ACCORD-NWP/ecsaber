/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Model;
  class ModelAuxIncrement;

// -----------------------------------------------------------------------------
/// ModelAuxControl class

class ModelAuxControl : public util::Printable,
                        private eckit::NonCopyable,
                        private util::ObjectCounter<ModelAuxControl> {
 public:
  static const std::string classname()
    {return "quench::ModelAuxControl";}

/// OOPS interface

// Constructors/destructor
  ModelAuxControl(const Geometry &,
                  const Model &,
                  const eckit::Configuration &)
    {}
  ModelAuxControl(const Geometry &,
                  const eckit::Configuration &)
    {}
  ModelAuxControl(const Geometry &,
                  const ModelAuxControl &)
    {}
  ModelAuxControl(const ModelAuxControl &,
                  const bool)
    {}
  ~ModelAuxControl()
    {}

// Basic operator
  ModelAuxControl & operator+=(const ModelAuxIncrement &)
    {return *this;}

// I/O and diagnostics
  void read(const eckit::Configuration &)
    {}
  void write(const eckit::Configuration &) const
    {}
  double norm() const
    {return 0.0;}

 private:
  void print(std::ostream & os) const
    {}
};

// -----------------------------------------------------------------------------

}  // namespace quench
