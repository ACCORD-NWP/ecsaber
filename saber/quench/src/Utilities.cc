/*
 * (C) Copyright 2025 Meteorologisk Institutt
 * 
 */

#include "src/Utilities.h"

#include <cmath>
#include <iomanip>

#include "util/dateFunctions.h"

namespace df = util::datefunctions;

namespace quench {

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

size_t dateTimeSerialSize(const util::DateTime & dateTime) {
  return 2;
}

// -----------------------------------------------------------------------------

void dateTimeSerialize(const util::DateTime & dateTime,
                       std::vector<double> & vect) {
  int year, month, day, hour, minute, second;
  dateTime.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  vect.push_back(static_cast<double>(df::dateToJulian(year, month, day)));
  vect.push_back(static_cast<double>(df::hmsToSeconds(hour, minute, second)));
}

// -----------------------------------------------------------------------------

void dateTimeDeserialize(util::DateTime & dateTime,
                         const std::vector<double> & vect,
                         size_t & current) {
  uint64_t date = std::lround(vect.at(current));
  int time = std::lround(vect.at(current+1));
  int year, month, day, hour, minute, second;
  df::julianToDate(date, year, month, day);
  df::secondToHms(time, hour, minute, second);
  dateTime = util::DateTime(year, month, day, hour, minute, second);
  current += 2;
}

// -----------------------------------------------------------------------------

}  // namespace quench
