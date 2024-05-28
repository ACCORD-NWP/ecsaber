/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "src/Fields.h"

#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class GeoVaLs;
  class Geometry;
  class Increment;
  class LinearModel;
  class Locations;
  class Model;
  class Variables;

// -----------------------------------------------------------------------------
/// State class

class State : public util::Printable,
              private util::ObjectCounter<State> {
 public:
  static const std::string classname()
    {return "quench::State";}

/// Constructor, destructor
  State(const Geometry &,
        const eckit::Configuration &);
  State(const Geometry & resol,
        const Model &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const LinearModel &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const State & other)
    : fields_(new Fields(*other.fields_, resol)) {}
  State(const Geometry & resol,
        const Model &,
        const State & other)
    : State(resol, other) {}
  State(const Geometry & resol,
        const Model &,
        const State & other,
        const eckit::Configuration &)
    : State(resol, other) {}
  State(const State & other)
    : fields_(new Fields(*other.fields_)) {}
  virtual ~State()
    {}

/// Assignment
  State & operator=(const State &);

/// Interpolate to observations locations
  void interpolate(const Locations & locs,
                   GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}

/// Force full fields
  void forceWith(const State & other,
                 const Variables & vars)
    {fields_->forceWith(*(other.fields_), vars);}

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const eckit::Configuration & config)
    {fields_->read(config);}
  void write(const eckit::Configuration & config) const
    {fields_->write(config);}
  double norm() const
    {return fields_->norm();}
  const util::DateTime & validTime() const
    {return fields_->time();}
  util::DateTime & validTime()
    {return fields_->time();}

/// Access to fields
  Fields & fields()
    {return *fields_;}
  const Fields & fields() const
    {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const
    {return fields_->geometry();}

/// ATLAS FieldSet accessor
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

/// Other
  void zero()
    {fields_->zero();}
  void accumul(const double & zz,
               const State & xx)
    {fields_->axpy(zz, xx.fields());}

 private:
  void print(std::ostream &) const;

  std::unique_ptr<Fields> fields_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
