/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/Variables.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {

// -----------------------------------------------------------------------------

JediVariables::JediVariables(const eckit::Configuration & conf, const std::string & name) {
  std::vector<std::string> vars;
  if (!conf.get(name, vars)) {
    Log::error() << name << " not found in " << conf << std::endl;
    throw eckit::BadParameter("Undefined variable: '" + name + "'");
  }
  vars_.reserve(vars.size());
  for (const auto & var : vars) {
    vars_.push_back(Variable(var));
  }
}

// -----------------------------------------------------------------------------

JediVariables::JediVariables(const std::vector<std::string> & vars) {
  vars_.reserve(vars.size());
  for (const auto & var : vars) {
    vars_.push_back(Variable(var));
  }
}

// -----------------------------------------------------------------------------

JediVariables::JediVariables(const std::vector<Variable> & vars) : vars_(vars) {
}

// -----------------------------------------------------------------------------

JediVariables::JediVariables(const eckit::Configuration & conf, const std::vector<std::string> & vars) {
  vars_.reserve(vars.size());
  for (const auto & var : vars) {
    if (conf.has(var)) {
      vars_.push_back(Variable(var, conf.getSubConfiguration(var)));
    } else {
      vars_.push_back(Variable(var));
    }
  }
}

// -----------------------------------------------------------------------------

const oops::Variable & JediVariables::operator[](const std::string & varname) const {
  const auto & it = std::find(vars_.begin(), vars_.end(), Variable(varname));
  if (it == vars_.end()) {
    Log::error() << "Could not find " << varname << " in variables list: " << *this << std::endl;
    throw eckit::BadParameter("Error accessing a non-existent variable", Here());
  }
  return *it;
}

// -----------------------------------------------------------------------------

oops::Variable & JediVariables::operator[](const std::string & varname) {
  const auto & it = std::find(vars_.begin(), vars_.end(), Variable(varname));
  if (it == vars_.end()) {
    Log::error() << "Could not find " << varname << " in variables list: " << *this << std::endl;
    throw eckit::BadParameter("Error accessing a non-existent variable", Here());
  }
  return *it;
}

// -----------------------------------------------------------------------------

const std::vector<std::string> JediVariables::variables() const {
  std::vector<std::string> vars;
  vars.reserve(vars_.size());
  for (const auto & var : vars_) {
    vars.push_back(var.name());
  }
  return vars;
}

// -----------------------------------------------------------------------------

bool JediVariables::has(const oops::Variable & var) const {
  return (std::find(vars_.begin(), vars_.end(), var) != vars_.end());
}

// -----------------------------------------------------------------------------

bool JediVariables::has(const std::string & varname) const {
  Variable var(varname);
  return this->has(var);
}

// -----------------------------------------------------------------------------

size_t JediVariables::find(const oops::Variable & var) const {
  const auto & it = std::find(vars_.begin(), vars_.end(), var);
  if (it == vars_.end()) {
    Log::error() << "Could not find " << var << " in variables list: " << *this << std::endl;
  }
  return (it - vars_.begin());
}

// -----------------------------------------------------------------------------

size_t JediVariables::find(const std::string & varname) const {
  Variable var(varname);
  return this->find(var);
}

// -----------------------------------------------------------------------------

void JediVariables::push_back(const Variable & var) {
  // Prevent adding the same variable twice
  if (auto it = std::find(vars_.begin(), vars_.end(), var); it == vars_.end()) {
    vars_.push_back(var);
  } else if ((it->getLevels() < 0) && (var.getLevels() >= 0)) {
    // This is required since negative levels are ignored in the '==' Variable comparison,
    // so we want to preserve the levels if they are set in the new variable.
    it->setLevels(var.getLevels());
  }
}

// -----------------------------------------------------------------------------

void JediVariables::push_back(const std::string & var) {
  push_back(Variable(var));
}

// -----------------------------------------------------------------------------

std::vector<Variable> JediVariables::asCanonical() const {
  std::vector<Variable> vars(vars_);
  std::sort(vars.begin(), vars.end());
  return vars;
}

// -----------------------------------------------------------------------------

JediVariables & JediVariables::operator+=(const JediVariables & rhs) {
  for (auto & var : rhs) {
    this->push_back(var);
  }
  return *this;
}

// -----------------------------------------------------------------------------

JediVariables & JediVariables::operator-=(const JediVariables & rhs) {
  for (auto & var : rhs.vars_) {
    vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  }
  return *this;
}

// -----------------------------------------------------------------------------

JediVariables & JediVariables::operator-=(const Variable & rhs) {
  vars_.erase(std::remove(vars_.begin(), vars_.end(), rhs), vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

bool JediVariables::operator==(const JediVariables & rhs) const {
  if (vars_.size() != rhs.vars_.size()) {
    return false;
  }
  std::vector<Variable> myvars = this->asCanonical();
  std::vector<Variable> othervars = rhs.asCanonical();
  for (size_t jvar = 0; jvar < myvars.size(); ++jvar) {
    if (myvars[jvar] != othervars[jvar]) {
      return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------

bool JediVariables::operator!=(const JediVariables & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

bool JediVariables::operator<=(const JediVariables & rhs) const {
  bool is_in_rhs = true;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    is_in_rhs = is_in_rhs && rhs.has(vars_[jj]);
  }
  return is_in_rhs;
}

// -----------------------------------------------------------------------------

void JediVariables::intersection(const JediVariables & rhs) {
  std::vector<Variable> myvars = this->asCanonical();
  std::vector<Variable> othervars = rhs.asCanonical();
  std::vector<Variable> commonvars;
  std::set_intersection(myvars.cbegin(), myvars.cend(),
                        othervars.cbegin(), othervars.cend(), std::back_inserter(commonvars));
  vars_ = commonvars;
}

// -----------------------------------------------------------------------------

void JediVariables::sort() {
  std::sort(vars_.begin(), vars_.end());
}

// -----------------------------------------------------------------------------

void JediVariables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
}

}  // namespace oops
