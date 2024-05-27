/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/serialisation/Stream.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/FieldSetOperations.h"

#include "quench/Interpolation.h"
#include "quench/Variables.h"

namespace quench {
  class Geometry;
  class GeoVaLs;
  class Locations;

// -----------------------------------------------------------------------------
/// Fields class

class Fields : public util::Printable,
               private util::ObjectCounter<Fields> {
 public:
  static const std::string classname()
    {return "quench::Fields";}

/// OOPS interface

// Constructors/destructor
  Fields(const Geometry &,
         const Variables &,
         const util::DateTime &);
  Fields(const Fields &,
         const Geometry &);
  Fields(const Fields &,
         const bool);
  Fields(const Fields &);
  ~Fields()
    {}

// Basic operators
  void zero();
  void constantValue(const double &);
  void constantValue(const eckit::Configuration &);
  Fields & operator=(const Fields &);
  Fields & operator+=(const Fields &);
  Fields & operator-=(const Fields &);
  Fields & operator*=(const double &);
  void axpy(const double &,
            const Fields &);
  double dot_product_with(const Fields &) const;
  void schur_product_with(const Fields &);
  void dirac(const eckit::Configuration &);
  void random();
  void diff(const Fields &,
            const Fields &);

// ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);
  const atlas::FieldSet & fieldSet() const
    {return fset_;}
  atlas::FieldSet & fieldSet()
    {return fset_;}

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const Geometry> geometry() const
    {return geom_;}
  const Variables & variables() const
    {return vars_;}
  const util::DateTime & time() const
    {return time_;}
  util::DateTime & time()
    {return time_;}
  void updateTime(const util::Duration & dt)
    {time_ += dt;}

/// Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &,
                   size_t &);

  atlas::FieldSet & fields() {return fset_;}

/// Grid interpolations
  static std::vector<quench::Interpolation>& interpolations();

 private:
  void print(std::ostream &) const;
  std::vector<quench::Interpolation>::iterator setupGridInterpolation(const Geometry &) const;

  std::shared_ptr<const Geometry> geom_;
  Variables vars_;
  util::DateTime time_;
  mutable atlas::FieldSet fset_;

/// ECSABER-specific interface
 public:
  double min(const Variables &) const;
  double max(const Variables &) const;
  void interpolate(const Locations &,
                   GeoVaLs &) const;
  void interpolateAD(const Locations &,
                     const GeoVaLs &);
  void forceWith(const Fields &,
                 const Variables &);
  void synchronizeFields();
  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Fields &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Fields &);

 private:
  std::vector<quench::Interpolation>::iterator setupObsInterpolation(const Locations &) const;
};

// -----------------------------------------------------------------------------

}  // namespace quench
