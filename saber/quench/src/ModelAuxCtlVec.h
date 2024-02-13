/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ModelAuxControl;
  class ModelAuxCovariance;

// -----------------------------------------------------------------------------
/// ModelAuxCtlVec Class

class ModelAuxCtlVec : public util::Printable {
 public:
  static const std::string classname()
    {return "quench::ModelAuxCtlVec";}

/// Constructor, destructor
  explicit ModelAuxCtlVec(const ModelAuxCovariance &)
    {}
  ModelAuxCtlVec(const ModelAuxCovariance &,
                 const ModelAuxCtlVec &)
    {}
  ModelAuxCtlVec(const ModelAuxCtlVec &,
                 const bool copy = true)
    {}
  ~ModelAuxCtlVec()
    {}

/// Linear algebra operators
  void zero()
    {}
  ModelAuxCtlVec & operator=(const ModelAuxCtlVec &)
    {return *this;}
  ModelAuxCtlVec & operator+=(const ModelAuxCtlVec &)
    {return *this;}
  ModelAuxCtlVec & operator-=(const ModelAuxCtlVec &)
    {return *this;}
  ModelAuxCtlVec & operator*=(const double &)
    {return *this;}
  void axpy(const double &,
            const ModelAuxCtlVec &)
    {}
  double dot_product_with(const ModelAuxCtlVec &) const
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
  void print(std::ostream &) const;

  const bool active_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace quench
