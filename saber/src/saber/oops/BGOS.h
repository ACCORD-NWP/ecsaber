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
#include "oops/base/Increment4D.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/Utilities.h"

namespace saber {

template <typename MODEL>
class BGOS : public oops::Application {
  using CostJo_ = oops::CostJo<MODEL>;
  using CtrlInc_ = oops::ControlIncrement<MODEL>;
  using CtrlVariable_ = oops::ControlVariable<MODEL>;
  using Departures_ = oops::Departures<MODEL>;
  using Dual_ = oops::DualVector<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using HMatrix_ = oops::HMatrix<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using State_ = oops::State<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  BGOS() {
    oops::instantiateCostFactory<MODEL>();
    saber::instantiateCovarFactory<MODEL>();
    oops::instantiateMinFactory<MODEL>();
    oops::instantiateObsErrorFactory<MODEL>();
    oops::instantiateTlmFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~BGOS() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    // Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    // Setup model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    // Setup cost function
    const eckit::LocalConfiguration cfConf(fullConfig, "cost_function");
    std::unique_ptr<oops::CostFunction<MODEL> > J(
        oops::CostFactory<MODEL>::create(cfConf, resol, model));
    oops::Log::trace() << "BGOS: cost function has been set up" << std::endl;

    // Initialize first guess from background
    const CtrlVariable_ xx(J->jb().getBackground());

    // Linearize cost function
    oops::PostProcessor<State_> post;
    std::vector<eckit::LocalConfiguration> iterconfs;
    fullConfig.get("variational.iteration", iterconfs);
    J->linearize(xx, iterconfs[0], post);

    // Get Jo
    const CostJo_ & Jo = dynamic_cast<const CostJo_ &>(J->jterm(0));

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

    // Setup linearized observation operator
    const HMatrix_ HMatrix(*J);

    // Setup dual vectors
    ASSERT(J->nterms() == 1);
    Dual_ mean;
    Dual_ var;
    mean.append(Jo.newDualVector());
    var.append(Jo.newDualVector());

    // Initialize mean and variance
    mean.zero();
    var.zero();

    // Convert ensemble of backgrounds perturbations from primal to dual space
    const bool nonlinear = fullConfig.getBool("use nonlinear observation operator", false);

    for (size_t ie = 0; ie < nens; ++ie) {
      // Get ensemble member
      Increment4D_ dx4D(J->jb().jbState());
      if (membersConfig.empty()) {
        // Randomize member as increment
        J->jb().jbState().randomize(dx4D);
      } else {
        // Read member as increment
        dx4D[dx4D.first()].read(membersConfig[ie]);
      }

      // Convert to observation space
      Dual_ dy;

      if (nonlinear) {
        // Use nonlinear observation operator

        // Re-evaluate J around new state
        CtrlVariable_ xxTmp(xx);
        if (!membersConfig.empty()) {
          xxTmp.state()[0].zero();       
        }
        xxTmp.state()[0] += dx4D[0];
        J->evaluate(xxTmp, eckit::LocalConfiguration());

        // Get departures
        dy.append(Jo.newDepartures());
      } else {
        // Use linearized observation operator

        // Setup control increment
        CtrlInc_ dx(J->jb());
        dx.state() = dx4D;

        // Allocate dual vector
        dy.append(Jo.newDualVector());

        // Apply H matrix to control increment
        HMatrix.multiply(dx, dy);
      }

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
