/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"

namespace util {

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration &,
               const size_t &);

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

// -----------------------------------------------------------------------------

void expandEnsembleTemplate(eckit::LocalConfiguration &,
                            const size_t &);

// -----------------------------------------------------------------------------

}  // namespace util
