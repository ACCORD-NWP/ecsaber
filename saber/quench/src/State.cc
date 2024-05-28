/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/State.h"

#include <vector>

#include "src/Fields.h"
#include "src/Increment.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

State::State(const Geometry & resol,
             const eckit::Configuration & file)
  : fields_() {
  oops::Log::trace() << classname() << "::State starting" << std::endl;

  const Variables vars(file);
  fields_.reset(new Fields(resol, vars, util::DateTime()));
  if (file.has("filepath")) {
    oops::Log::info() << "Info     : Create state from file" << std::endl;
    fields_->read(file);
  } else {
    oops::Log::info() << "Info     : Create empty state" << std::endl;
    if (file.has("constant value")) {
      fields_->constantValue(file.getDouble("constant value"));
    } else if (file.has("constant group-specific value")) {
      fields_->constantValue(file);
    } else {
      fields_->zero();
    }
  }
  const util::DateTime vt(file.getString("date"));
  fields_->time() = vt;

  oops::Log::trace() << classname() << "::State done" << std::endl;
}

// -----------------------------------------------------------------------------

State & State::operator=(const State & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;

  oops::Log::trace() << classname() << "::operator= done" << std::endl;
}

// -----------------------------------------------------------------------------

State & State::operator+=(const Increment & dx) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  *fields_+=dx.fields();

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void State::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << std::endl << "Valid time:" << this->validTime();
  os << *fields_;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
