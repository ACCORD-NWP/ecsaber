/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/Fields.h"
#include "src/StateECInc.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Increment;

// -----------------------------------------------------------------------------
/// State class

class State : public util::Printable,
              private util::ObjectCounter<State> {
 public:
  static const std::string classname()
    {return "quench::State";}

  // Constructors
  State(const Geometry &,
        const Variables &,
        const util::DateTime &);
  State(const Geometry &,
        const eckit::Configuration &);
  State(const Geometry & resol,
        const State & other)
    : fields_(new Fields(*other.fields_, resol)) {}
  State(const Variables & vars,
        const State & other)
    : fields_(new Fields(*other.fields_)) {}
  State(const State & other)
    : fields_(new Fields(*other.fields_)) {}

  // Assignment
  State & operator=(const State &);

  // Interactions with Increment
  State & operator+=(const Increment &);

  // I/O and diagnostics
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

  // Access to fields
  Fields & fields()
    {return *fields_;}
  const Fields & fields() const
    {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const
    {return fields_->geometry();}

  // ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet & fset) const
    {fields_->toFieldSet(fset);}
  void fromFieldSet(const atlas::FieldSet & fset)
    {fields_->fromFieldSet(fset);}

  // Other
  void zero()
    {fields_->zero();}
  void accumul(const double & zz,
               const State & xx)
    {fields_->axpy(zz, xx.fields());}
  const Variables & variables() const
    {return fields_->variables();}

  // Serialization
  size_t serialSize() const
    {return fields_->serialSize();}
  void serialize(std::vector<double> & vect) const
    {fields_->serialize(vect);}
  void deserialize(const std::vector<double> & vect,
                   size_t & index)
    {fields_->deserialize(vect, index);}

 private:
  // Print
  void print(std::ostream &) const;

  // Fields
  std::unique_ptr<Fields> fields_;

/// ECSABER-specific definitions
#include "src/StateECDef.h"
};

// -----------------------------------------------------------------------------

}  // namespace quench
