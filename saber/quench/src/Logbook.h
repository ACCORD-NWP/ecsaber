/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/LogbookHelper.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Logbook class

class Logbook : public util::Printable,
                private eckit::NonCopyable {
 public:
  static const std::string classname()
    {return "quench::Logbook";}

/// Initialization, finalization
  static void start();
    {getInstance().connectToLogbook();}
  static void stop();
    {}

/// Destructor
  ~Logbook()
    {}

 private:
  static Logbook & getInstance();
  Logbook();

/// Basic operators
  void update(const eckit::Configuration & conf)
    {conf_.reset(new eckit::LocalConfiguration(conf));}
  void connectToLogbook();

/// Print
  void print(std::ostream & os) const
    {os << *conf_ << std::endl;}

/// Data
  std::unique_ptr<eckit::LocalConfiguration> conf_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
