/*
 * (C) Copyright 2025 Meteorologisk Institutt
 * 
 */

#pragma once

#include <vector>

#include "oops/util/DateTime.h"

namespace quench {

// -----------------------------------------------------------------------------

std::string dateTimeToStringIO(const util::DateTime &);

// -----------------------------------------------------------------------------

size_t dateTimeSerialSize(const util::DateTime &);

// -----------------------------------------------------------------------------

void dateTimeSerialize(const util::DateTime &,
                       std::vector<double> &);

// -----------------------------------------------------------------------------

void dateTimeDeserialize(util::DateTime &,
                         const std::vector<double> &,
                         size_t &);

// -----------------------------------------------------------------------------

}  // namespace quench
