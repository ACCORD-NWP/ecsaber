/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/memory/NonCopyable.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Geometry.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ModelAuxControl;
  class Fields;
  class State;

// -----------------------------------------------------------------------------
///  Model class

class Model: public util::Printable,
             private eckit::NonCopyable,
             private util::ObjectCounter<Model> {
 public:
  static const std::string classname()
    {return "quench::Model";}

  Model(const Geometry &,
        const eckit::Configuration &)
    {}
  Model(const Model &)
    {}
  ~Model()
    {}

/// Prepare model integration
  void initialize(State &,
                  const ModelAuxControl &) const
    {}

/// Model integration
  void step(State &,
            const ModelAuxControl &) const
    {}
  int saveTrajectory(State &,
                     const ModelAuxControl &) const
    {return 0;}

/// Finish model integration
  void finalize(State &) const
    {}

/// Utilities
  const util::Duration & timeResolution() const
    {return timeResolution_;}

 private:
  void print(std::ostream &) const
    {}

  const util::Duration timeResolution_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
