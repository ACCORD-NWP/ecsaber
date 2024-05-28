/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ObsAuxControl.h"

#include <iostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "util/Logger.h"

#include "src/ObsAuxIncrement.h"


namespace quench {

// -----------------------------------------------------------------------------

static oops::ObsAuxControlMaker<Traits, ObsAuxControl> makesObsAuxControlDefault_("default");

// -----------------------------------------------------------------------------

ObsAuxControl::ObsAuxControl(const ObsSpace &,
                             const eckit::Configuration & conf)
  : bias_(0.0), active_(false), obsSpaceKey_(0), key_(0) {
  oops::Log::trace() << classname() << "::ObsAuxControl starting" << std::endl;

  if (conf.has("bias")) {
    bias_ = conf.getDouble("bias");
    active_ = true;
  }

  oops::Log::trace() << classname() << "::ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAuxControl::ObsAuxControl(const ObsAuxControl & other,
                             const bool copy)
  : bias_(0.0), active_(other.active_), obsSpaceKey_(other.obsSpaceKey_), key_(other.key_) {
  oops::Log::trace() << classname() << "::ObsAuxControl starting" << std::endl;

  if (active_ && copy) {
    bias_ = other.bias_;
  }

  oops::Log::trace() << classname() << "::ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsAuxControlBase<Traits> & ObsAuxControl::operator+=(
  const oops::ObsAuxIncrementBase<Traits> & dx) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  if (active_) {
    const ObsAuxIncrement * pdx = dynamic_cast<const ObsAuxIncrement*>(&dx);
    bias_ += (*pdx).value();
  }

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsAuxControl::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  if (active_) {
    os << std::endl << "ObsAuxControl = " << bias_;
  }

  oops::Log::trace() << classname() << "::print dp,e" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
