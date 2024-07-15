/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class Configuration;
}

namespace quench {
  class Increment;

// -----------------------------------------------------------------------------
/// HorizScaleDecomposition class

class HorizScaleDecomposition: public util::Printable,
                               private eckit::NonCopyable,
                               private util::ObjectCounter<HorizScaleDecomposition> {
 public:
  static const std::string classname()
    {return "quench::HorizScaleDecomposition";}

/// OOPS interface
  explicit HorizScaleDecomposition(const eckit::Configuration &)
    {throw eckit::NotImplemented(Here());}

  void collect(const Increment &,
               const int &) const
    {throw eckit::NotImplemented(Here());}
  void decompose(const Increment &,
                 const int &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream & os) const
    {throw eckit::NotImplemented(Here());}
};
// -----------------------------------------------------------------------------

}  // namespace quench
