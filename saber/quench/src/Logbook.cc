/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Logbook.h"

namespace quench {

// -----------------------------------------------------------------------------

Logbook & Logbook::getInstance() {
  oops::Log::trace() << "Logbook::getInstance starting" << std::endl;

  static Logbook theLogbook;

  oops::Log::trace() << "Logbook::getInstance done" << std::endl;
  return theLogbook;
}

// -----------------------------------------------------------------------------

Logbook::Logbook()
  : conf_(new eckit::LocalConfiguration()) {
  oops::Log::trace() << "Logbook::Logbook done" << std::endl;
}

// -----------------------------------------------------------------------------

void Logbook::connectToLogbook() {
  oops::Log::trace() << "Logbook::connectToLogbook done" << std::endl;
  util::LogbookHelper::connectCallback([this]
                                       (const eckit::Configuration & conf)
                                       {return update(conf);});
}

// -----------------------------------------------------------------------------

}  // namespace quench
