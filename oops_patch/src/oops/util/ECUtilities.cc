/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "oops/util/ECUtilities.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <string>

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/ConfigFunctions.h"

namespace util {

// -----------------------------------------------------------------------------
// Timestamp
// -----------------------------------------------------------------------------

static std::chrono::steady_clock::time_point start_time(std::chrono::steady_clock::now());

double timeStamp() {
  const std::chrono::steady_clock::time_point t(std::chrono::steady_clock::now());
  return std::chrono::duration<double>(t - start_time).count();
}

// -----------------------------------------------------------------------------
// DateTime
// -----------------------------------------------------------------------------

std::string dateTimeToStringIO(const util::DateTime & dateTime) {
  int year, month, day, hour, minute, second;
  dateTime.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);

  std::ostringstream os;
  os << std::setfill('0');
  os << std::setw(4) << year;
  os << std::setw(2) << month;
  os << std::setw(2) << day;
  os.put('T');
  os << std::setw(2) << hour;
  os << std::setw(2) << minute;
  os << std::setw(2) << second;
  os.put('Z');
  return os.str();
}

// -----------------------------------------------------------------------------

}  // namespace util
