/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class ModelAuxControl;
  class ModelAuxCovariance;

// -----------------------------------------------------------------------------
/// ModelAuxIncrement class

class ModelAuxIncrement : public util::Printable {
 public:
  static const std::string classname()
    {return "quench::ModelAuxIncrement";}

/// OOPS interface

// Constructors/destructor
  ModelAuxIncrement(const Geometry &,
                    const eckit::Configuration &)
    {}
  ModelAuxIncrement(const ModelAuxIncrement &,
                    const bool = true)
    {}
  ModelAuxIncrement(const ModelAuxIncrement &,
                    const eckit::Configuration &)
    {}
  ~ModelAuxIncrement()
    {}

// Linear algebra operators
  void diff(const ModelAuxControl &,
            const ModelAuxControl &)
    {}
  void zero()
    {}
  ModelAuxIncrement & operator=(const ModelAuxIncrement &)
    {return *this;}
  ModelAuxIncrement & operator+=(const ModelAuxIncrement &)
    {return *this;}
  ModelAuxIncrement & operator-=(const ModelAuxIncrement &)
    {return *this;}
  ModelAuxIncrement & operator*=(const double &)
    {return *this;}
  void axpy(const double,
            const ModelAuxIncrement &)
    {}
  double dot_product_with(const ModelAuxIncrement &) const
    {return 0.0;}

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
