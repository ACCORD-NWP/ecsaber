/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const Variables & vars,
                     const util::DateTime &,
                     const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)) {
  oops::Log::trace() << classname() << "::Increment starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << classname() << "::Increment done" << std::endl;
}


// -----------------------------------------------------------------------------

eckit::Stream & operator<<(eckit::Stream & s,
                           const Increment & dx) {
  oops::Log::trace() << "Increment::operator<< starting" << std::endl;

  s << dx.fields();

  oops::Log::trace() << "Increment::operator<< done" << std::endl;
  return s;
}

// -----------------------------------------------------------------------------

eckit::Stream & operator>>(eckit::Stream & s,
                           Increment & dx) {
  oops::Log::trace() << "Increment::operator>> starting" << std::endl;

  s >> dx.fields();

  oops::Log::trace() << "Increment::operator>> done" << std::endl;
  return s;
}

// -----------------------------------------------------------------------------
