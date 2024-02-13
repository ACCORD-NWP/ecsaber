/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Increment.h"

#include <vector>

#include "atlas/field.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const Variables & vars,
                     const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)) {
  oops::Log::trace() << "Increment::Increment starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << "Increment::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const Variables & vars,
                     const util::DateTime &,
                     const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)) {
  oops::Log::trace() << "Increment::Increment starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << "Increment::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const Increment & other)
  : fields_() {
  oops::Log::trace() << "Increment::Increment starting" << std::endl;

  throw eckit::NotImplemented(Here());

  oops::Log::trace() << "Increment::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment::Increment(const Increment & other,
                     const bool copy)
  : fields_(new Fields(*other.fields_, copy)) {
  oops::Log::trace() << "Increment::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::diff(const State & x1,
                     const State & x2) {
  oops::Log::trace() << "Increment::diff starting" << std::endl;

  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());

  oops::Log::trace() << "Increment::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator=(const Increment & rhs) {
  oops::Log::trace() << "Increment::operator= starting" << std::endl;

  *fields_ = *rhs.fields_;

  oops::Log::trace() << "Increment::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator+=(const Increment & dx) {
  oops::Log::trace() << "Increment::operator+= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;

  oops::Log::trace() << "Increment::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator-=(const Increment & dx) {
  oops::Log::trace() << "Increment::operator-= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;

  oops::Log::trace() << "Increment::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator*=(const double & zz) {
  oops::Log::trace() << "Increment::operator*= starting" << std::endl;

  *fields_ *= zz;

  oops::Log::trace() << "Increment::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void Increment::zero(const util::DateTime & vt) {
  oops::Log::trace() << "Increment::zero starting" << std::endl;

  fields_->zero();
  fields_->time() = vt;

  oops::Log::trace() << "Increment::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::axpy(const double & zz,
                     const Increment & dx,
                     const bool check) {
  oops::Log::trace() << "Increment::axpy starting" << std::endl;

  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);

  oops::Log::trace() << "Increment::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::print(std::ostream & os) const {
  oops::Log::trace() << "Increment::print starting" << std::endl;

  os << std::endl << "Valid time:" << this->validTime();
  os << *fields_;

  oops::Log::trace() << "Increment::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
