/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/serialisation/Stream.h"

#include "oops/base/GeneralizedDepartures.h"

#include "quench/Fields.h"
#include "quench/State.h"

#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class Duration;
}

namespace quench {
  class Geometry;
  class GeoVaLs;
  class Locations;
  class ModelAuxIncrement;
  class Variables;

// -----------------------------------------------------------------------------
/// Increment class

class Increment : public oops::GeneralizedDepartures,
                  public util::Printable,
                  private util::ObjectCounter<Increment> {
 public:
  static const std::string classname()
    {return "quench::Increment";}

/// OOPS interface

// Constructors/destructor
  Increment(const Geometry &,
            const Variables &,
            const util::DateTime &);
  Increment(const Geometry &,
            const Variables &,
            const util::DateTime &,
            const util::DateTime &);
  Increment(const Geometry &,
            const Increment &);
  Increment(const Increment &,
            const bool);
  Increment(const Increment &)
    {}
  virtual ~Increment()
    {}

// Basic operators
  void diff(const State &,
            const State &);
  void zero()
    {fields_->zero();}
  void zero(const util::DateTime &);
  void dirac(const eckit::Configuration & config)
    {fields_->dirac(config);}
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &,
            const Increment &,
            const bool check = true);
  double dot_product_with(const Increment & dx) const
    {return fields_->dot_product_with(*dx.fields_);}
  void schur_product_with(const Increment & dx)
    {fields_->schur_product_with(*dx.fields_);}
  void random()
    {fields_->random();}

// Interpolate to observations locations
  void interpolateTL(const Locations & locs,
                     GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}
  void interpolateAD(const Locations & locs,
                     const GeoVaLs & gv)
    {fields_->interpolateAD(locs, gv);}

// I/O and diagnostics
  void read(const eckit::Configuration & config)
    {fields_->read(config);}
  void write(const eckit::Configuration & config) const
    {fields_->write(config);}
  double norm() const
    {return fields_->norm();}
  double max(const Variables & var) const
    {return fields_->max(var);}
  double min(const Variables & var) const
    {return fields_->min(var);}
  const util::DateTime & validTime() const
    {return fields_->time();}
  util::DateTime & validTime()
    {return fields_->time();}
  void updateTime(const util::Duration & dt)
    {fields_->time() += dt;}

/// Access to fields
  Fields & fields()
    {return *fields_;}
  const Fields & fields() const
    {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const
    {return fields_->geometry();}

/// Other
  void accumul(const double & zz,
               const State & xx)
    {fields_->axpy(zz, xx.fields());}


/// Serialization
  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Increment &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Increment &);

/// ECSABER interface

// ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet & fset) const
    {fields_->toFieldSet(fset);}
  void fromFieldSet(const atlas::FieldSet & fset)
    {fields_->fromFieldSet(fset);}
  const atlas::FieldSet & fieldSet() const
    {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet()
    {return fields_->fieldSet();}
  void synchronizeFields()
    {fields_->synchronizeFields();}

 private:
  void print(std::ostream &) const;

  std::unique_ptr<Fields> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
