/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace quench {
  class Increment;
  class ObsSpace;
  class State;
  class Variables;

// -----------------------------------------------------------------------------
/// GeoVaLs class

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname()
    {return "quench::GeoVaLs";}

/// OOPS interface

// Constructors/destructor
  GeoVaLs(const ObsSpace &,
          const Variables &,
          const Increment &,
          const util::DateTime &,
          const util::DateTime &);
  GeoVaLs(const ObsSpace &,
          const Variables &,
          const State &,
          const util::DateTime &,
          const util::DateTime &);
  ~GeoVaLs()
    {}

// Basic operators
  void zero();
  void random();
  double dot_product_with(const GeoVaLs &) const;

// Read/write
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

/// Local
  size_t obsIndex(const size_t & io) const
    {return iobs_[io];}
  const atlas::FieldSet & fieldSet() const
    {return gvFieldSet_;}
  atlas::FieldSet & fieldSet()
    {return gvFieldSet_;}

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  const ObsSpace & obsSpace_;
  std::vector<size_t> iobs_;
  atlas::FieldSet gvFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
