/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

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

  GeoVaLs(const ObsSpace &,
          const Variables &,
          const Increment &,
          const util::DateTime &,
          const util::DateTime &)
    {throw eckit::NotImplemented(Here());}
  GeoVaLs(const ObsSpace &,
          const Variables &,
          const State &,
          const util::DateTime &,
          const util::DateTime &)
    {throw eckit::NotImplemented(Here());}

  void zero()
    {throw eckit::NotImplemented(Here());}
  void random()
    {throw eckit::NotImplemented(Here());}
  double dot_product_with(const GeoVaLs &) const
    {throw eckit::NotImplemented(Here()); return 0.0;}
  void read(const eckit::Configuration &)
    {throw eckit::NotImplemented(Here());}
  void write(const eckit::Configuration &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream &) const;
};

// -----------------------------------------------------------------------------

}  // namespace quench
