/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "util/Printable.h"

#include "src/Geometry.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ModelAuxControl;
  class ModelAuxCovariance;
  class Geometry;

// -----------------------------------------------------------------------------
/// ModelAuxIncrement class

class ModelAuxIncrement : public util::Printable {
 public:
  static const std::string classname()
    {return "quench::ModelAuxIncrement";}

/// Constructor, destructor
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

/// Linear algebra operators
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

/// I/O and diagnostics
  void read(const eckit::Configuration &)
    {}
  void write(const eckit::Configuration &) const
    {}
  double norm() const
    {return 0.0;}

/// Status
  const bool & isActive() const
    {return active_;}

 private:
  void print(std::ostream & os) const
    {}

  const bool active_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace quench
