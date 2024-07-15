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
#include "oops/assimilation/ControlObsVector.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/Departures.h"
#include "oops/base/Increment4D.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/Observations.h"
#include "oops/base/ObservationSpaces.h"
#include "oops/base/Observer.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"

namespace saber {

template <typename MODEL> 
class BGOS : public oops::Application {
  using ControlObsVector_ = oops::ControlObsVector<MODEL>;
  using CtrlInc_ = oops::ControlIncrement<MODEL>;
  using CtrlVariable_ = oops::ControlVariable<MODEL>;
  using Departures_ = oops::Departures<MODEL>;
  using Dual_ = oops::DualVector<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using HMatrix_ = oops::HMatrix<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using ModelAux_ = oops::ModelAuxControl<MODEL>;
  using ObsAuxCtrls_ = oops::ObsAuxControls<MODEL>;
  using Observations_ = oops::Observations<MODEL>;
  using ObsSpaces_ = oops::ObservationSpaces<MODEL>;
  using ObsOperators_ = oops::ObsOperators<MODEL>;
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

    // Define linearizing state
    const State_ xb(xx.state()[0]);

    // Window times
    const eckit::LocalConfiguration obsconf(fullConfig, "cost_function.Jo");
    const util::DateTime winbgn(fullConfig.getString("cost_function.window_begin"));
    const util::Duration winlen(fullConfig.getString("cost_function.window_length"));

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

    // Setup augmented state
    const ModelAux_ moderr(resol, model, eckit::LocalConfiguration());

    // Get observations space
    const ObsSpaces_ & obsdb = xx.obsVar().obspaces();

    // Setup observations
    const ObsOperators_ hop(obsdb);

    // Setup observation bias
    const ObsAuxCtrls_ ybias(obsdb, obsconf);

    // Setup linearized observation operator
    const HMatrix_ HMatrix(*J);

    // Setup dual vectors
    ASSERT(J->nterms() == 1);
    Dual_ mean;
    Dual_ var;
    mean.append(J->jterm(0).newDualVector());
    var.append(J->jterm(0).newDualVector());

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
      std::vector<std::shared_ptr<ControlObsVector_>> controlObsVecSharedVec;

      if (nonlinear) {
        // Use nonlinear observation operator

        // Add increment to state
        State_ xTmp(xb);
        if (!membersConfig.empty()) {
          xTmp.zero();       
        }
        xTmp += dx4D[0];

        // Setup forecast outputs
        oops::PostProcessor<State_> postH;
        std::shared_ptr<oops::Observer<MODEL, State_> > pobs(
          new oops::Observer<MODEL, State_>(obsdb, hop, ybias));
        postH.enrollProcessor(pobs);

        // Run forecast
        model.forecast(xTmp, moderr, winlen, postH);

        // Get observations
        std::unique_ptr<Observations_> yobs(pobs->release());

        // Move observation data to dual vector
        for (size_t jj = 0; jj < obsdb.size(); ++jj) {
          std::shared_ptr<ControlObsVector_> controlObsVec(new ControlObsVector_((*yobs)[jj]));
          controlObsVecSharedVec.push_back(controlObsVec);
        }
        dy.append(new Departures_(controlObsVecSharedVec));
      } else {
        // Use linearized observation operator

        // Setup control increment
        CtrlInc_ dx(J->jb());
        dx.state() = dx4D;

        // Allocate dual vector
        dy.append(J->jterm(0).newDualVector());

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
