/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace quench {

// -----------------------------------------------------------------------------
///  Variables class

class Variables : public oops::JediVariables,
                  private util::ObjectCounter<Variables> {
  using quenchVariables = oops::JediVariables;
  using quenchVariables::quenchVariables;

 public:
  static const std::string classname()
    {return "quench::Variables";}

// Extra constructor
  explicit Variables(const eckit::Configuration &);

// Extra accessor
  std::vector<std::string> variablesList() const
    {return this->variables();}
};

// -----------------------------------------------------------------------------

}  // namespace quench
