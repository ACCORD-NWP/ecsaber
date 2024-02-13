/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/interface/GenericMatrix.h"

#include "src/Traits.h"

#include "util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Hybrid weight class

class HybridWeight : public oops::GenericMatrix<Traits>,
                     private util::ObjectCounter<HybridWeight> {
 public:
  static const std::string classname()
    {return "quench::HybridWeight";}

  explicit HybridWeight(const eckit::Configuration &);
  ~HybridWeight()
    {}

  void multiply(State &) const
    {throw eckit::NotImplemented(Here());}
  void multiply(Increment & dx) const
    {dx *= wgt_;}
  void inverseMultiply(Increment &) const
    {throw eckit::NotImplemented(Here());}
  void transposeInverseMultiply(Increment &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream & os) const
    {os << std::endl << "Hybrid weight:" << wgt_;}

  double wgt_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
