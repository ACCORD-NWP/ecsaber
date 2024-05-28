/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/LinearModelBase.h"

#include "src/Traits.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {

// -----------------------------------------------------------------------------
/// LinearModel class

class LinearModel: public oops::LinearModelBase<Traits>,
           private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname()
    {return "quench::LinearModel";}

/// OOPS interface

// Constructors/destructor
  LinearModel(const Geometry &,
              const eckit::Configuration &);
  ~LinearModel()
    {}

// Model trajectory computation
  void setTrajectory(State &,
                     const Model &,
                     const ModelAuxControl &)
    {}

// Run TLM and its adjoint
  void initializeTL(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepTL(Increment &,
              const ModelAuxIncrement &) const;
  void finalizeTL(Increment &) const
    {}

  void initializeAD(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepAD(Increment &,
              ModelAuxIncrement &) const;
  void finalizeAD(Increment &) const
    {}

// Other utilities
  const util::Duration & timeResolution() const
    {return timeResolution_;}

 private:
  void print(std::ostream &) const
    {throw eckit::NotImplemented(Here());}

  const util::Duration timeResolution_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
