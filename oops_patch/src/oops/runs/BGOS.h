/*
 * (C) Copyright 2024 Norwegian Meteorological Institute.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/runs/Application.h"
#include "oops/runs/BGOS_impl.h"

namespace oops {

template <typename MODEL>
class BGOS : public Application {
 public:
  // -----------------------------------------------------------------------------
  BGOS() {
    instantiateCostFactory<MODEL>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL>();
    instantiateObsErrorFactory<MODEL>();
    instantiateTlmFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~BGOS() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    // Execute BGOS implementation
    BGOS_impl<MODEL>(fullConfig);

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::BGOS<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops
