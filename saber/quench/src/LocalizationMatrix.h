/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class Geometry;
  class Increment;
  class IncrEnsCtlVec;
  class Variables;

// -----------------------------------------------------------------------------
/// LocalizationMatrix class

class LocalizationMatrix: public util::Printable,
                          private eckit::NonCopyable,
                          private util::ObjectCounter<LocalizationMatrix> {
 public:
  static const std::string classname()
    {return "quench::LocalizationMatrix";}

  LocalizationMatrix(const Geometry &,
                     const Variables &,
                     const eckit::Configuration &)
    {throw eckit::NotImplemented(Here());}
  ~LocalizationMatrix()
    {}

  void multiply(Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrt(const IncrEnsCtlVec &,
                    Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrtTrans(const Increment &,
                         IncrEnsCtlVec &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream &) const
    {throw eckit::NotImplemented(Here());}
};
// -----------------------------------------------------------------------------

}  // namespace quench
