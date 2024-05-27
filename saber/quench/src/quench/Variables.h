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
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {

// -----------------------------------------------------------------------------
///  Variables class

class Variables : public util::Printable,
                  private util::ObjectCounter<Variables> {
 public:
  static const std::string classname()
    {return "quench::Variables";}

/// OOPS interface

  explicit Variables(const eckit::Configuration &);
  Variables(const Variables &);
  ~Variables() {}

  std::vector<std::string> variablesList() const
    {return vars_.variables();}
  size_t size() const
    {return vars_.size();}

/// Local

  explicit Variables(const std::vector<std::string> &);
  std::vector<std::string> variables() const
    {return vars_.variables();}
  const oops::JediVariables & toJediVariables() const
    {return vars_;}

 private:
  void print(std::ostream &) const;

  oops::JediVariables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
