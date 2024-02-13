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
  class Gom;
  class Geometry;
  class Increment;
  class Locations;
  class Model;
  class Tlm;
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
  State(const Geometry &,
        const Model &,
        const eckit::Configuration &);
  State(const Geometry &,
        const Tlm &,
        const eckit::Configuration &);
  State(const Geometry &,
        const State &);
  State(const Geometry &,
        const Model &,
        const State &);
  State(const Geometry &,
        const Model &,
        const State &,
        const eckit::Configuration &);
  State(const State &);
  virtual ~State()
    {}

/// Assignment
  State & operator=(const State &);

/// Interpolate to observation location
  void interpolate(const Locations &,
                   Gom &) const
    {throw eckit::NotImplemented(Here());}
/// Interpolate full fields
  void forceWith(const State &,
                 const Variables &)
    {throw eckit::NotImplemented(Here());}

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
  void activateModel()
    {throw eckit::NotImplemented(Here());}
  void deactivateModel()
    {throw eckit::NotImplemented(Here());}
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
