/*
 * (C) Copyright 2024 Norwegian Meteorological Institute.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <fstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/Departures.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"

namespace saber {

template <typename MODEL>
class BGOS : public oops::Application {
  using CtrlInc_ = oops::ControlIncrement<MODEL>;
  using CtrlVariable_ = oops::ControlVariable<MODEL>;
  using Departures_ = oops::Departures<MODEL>;
  using Dual_ = oops::DualVector<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using HMatrix_ = oops::HMatrix<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using State_ = oops::State<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  BGOS() {
    oops::instantiateCostFactory<MODEL>();
    oops::instantiateCovarFactory<MODEL>();
    oops::instantiateMinFactory<MODEL>();
    oops::instantiateObsErrorFactory<MODEL>();
    oops::instantiateTlmFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~BGOS() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    //  Setup model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    // Setup cost function
    const eckit::LocalConfiguration cfConf(fullConfig, "cost_function");
    std::unique_ptr<oops::CostFunction<MODEL> > J(
        oops::CostFactory<MODEL>::create(cfConf, resol, model));
    oops::Log::trace() << "BGOS: cost function has been set up" << std::endl;

    // Initialize first guess from background
    CtrlVariable_ xx(J->jb().getBackground());

    // Linearize cost function
    oops::PostProcessor<State_> post;
    std::vector<eckit::LocalConfiguration> iterconfs;
    fullConfig.get("variational.iteration", iterconfs);
    J->linearize(xx, iterconfs[0], post);

    // Define linearizing state
    const State_ xb(xx.state()[0]);

    // Setup H matrix
    HMatrix_ HMatrix(*J);

    // Convert ensemble of backgrounds perturbations from primal to dual space

    // Setup ControlIncrement and DualVector
    CtrlInc_ dx(J->jb());
    Dual_ dy;
    for (unsigned jj = 0; jj < J->nterms(); ++jj) {
      dy.append(J->jterm(jj).newDualVector());
    }

    // Setup mean and variance
    Dual_ mean(dy);
    Dual_ var(dy);
    mean.zero();
    var.zero();

    // Setup ensemble configuration (if present) or randomization size
    std::vector<eckit::LocalConfiguration> membersConfig;
    size_t nens;
    if (fullConfig.has("ensemble")) {
      eckit::LocalConfiguration membersConfigTemplate;
      fullConfig.get("ensemble", membersConfigTemplate);
      nens = membersConfigTemplate.getInt("members");
      expandEnsembleTemplate(membersConfigTemplate, nens);
      membersConfig = membersConfigTemplate.getSubConfigurations("state");
    } else {
      fullConfig.get("randomization size", nens);
    }
    ASSERT(nens > 1);

    for (size_t ie = 0; ie < nens; ++ie) {
      if (membersConfig.empty()) {
        // Randomize member as increment
        J->jb().jbState().randomize(dx.state());
      } else {
        // Read member as increment
        dx.state()[dx.state().first()].read(membersConfig[ie]);
      }

      // Apply H matrix to increment
      HMatrix.multiply(dx, dy);

      // Remove mean
      dy -= mean;

      // Update variance
      if (ie > 0) {
        const double fac = static_cast<double>(ie)/static_cast<double>(ie+1);
        Dual_ dy2(dy);
        dy2.schur_product_with(dy);
        var.axpy(fac, dy2);
      }

      // Update mean
      const double fac = 1.0/static_cast<double>(ie+1);
      mean.axpy(fac, dy);
    }

    // Normalize variance
    const double fac = 1.0/static_cast<double>(nens-1);
    var *= fac;

    // Write variance
    const eckit::LocalConfiguration varianceConfig(fullConfig, "variance in observation space");
    var.write(varianceConfig);

    // Print global norm
    oops::Log::test() << "Variance in observation space norm: " << 
      std::sqrt(var.dot_product_with(var)) << std::endl;

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "saber::BGOS<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace saber