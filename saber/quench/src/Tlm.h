/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/LinearModelBase.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Traits.h"

namespace eckit {
  class Configuration;
}

namespace quench {

// -----------------------------------------------------------------------------
/// Tlm class

class Tlm: public oops::LinearModelBase<Traits>,
           private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname()
    {return "quench::Tlm";}

  Tlm(const Geometry &,
      const eckit::Configuration &)
    {}
  ~Tlm()
    {}

/// Model trajectory computation
  void setTrajectory(State &,
                     const Model &,
                     const ModelAuxControl &)
    {}

/// Run TLM and its adjoint
  void initializeTL(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepTL(Increment &,
              const ModelAuxIncrement &) const
    {}
  void finalizeTL(Increment &) const
    {}

  void initializeAD(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepAD(Increment &,
              ModelAuxIncrement &) const
    {}
  void finalizeAD(Increment &) const
    {}

/// Other utilities
  const util::Duration & timeResolution() const
    {}

 private:
  void print(std::ostream &) const
    {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
