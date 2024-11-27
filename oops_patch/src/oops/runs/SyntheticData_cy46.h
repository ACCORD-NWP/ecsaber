/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/Departures.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/Observations.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ECUtilities.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SyntheticData : public Application {
  using CostJo_ = CostJo<MODEL>;
  using CovarianceBase_ = ModelSpaceCovarianceBase<MODEL>;
  using CovarianceFactory_ = CovarianceFactory<MODEL>;
  using CtrlVariable_ = ControlVariable<MODEL>;
  using Departures_ = Departures<MODEL>;
  using Geometry_ = Geometry<MODEL>;
  using HMatrix_ = HMatrix<MODEL>;
  using Increment_ = Increment<MODEL>;
  using Model_ = Model<MODEL>;
  using Observations_ = Observations<MODEL>;
  using State_ = State<MODEL>;
  using Variables_ = Variables<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  SyntheticData() {
    instantiateCostFactory<MODEL>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL>();
    instantiateObsErrorFactory<MODEL>();
    instantiateTlmFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~SyntheticData() {}
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
    std::unique_ptr<CostFunction<MODEL> > J(
        CostFactory<MODEL>::create(cfConf, resol, model));
    Log::trace() << "SyntheticData: cost function has been set up" << std::endl;

    // Initialize first guess from background
    const CtrlVariable_ xx(J->jb().getBackground());

    // Linearize cost function
    PostProcessor<State_> post;
    std::vector<eckit::LocalConfiguration> iterconfs;
    fullConfig.get("variational.iteration", iterconfs);
    J->linearize(xx, iterconfs[0], post);

    // Get observations equivalent from the truth
    const CostJo_ & Jo = dynamic_cast<const CostJo_ &>(J->jterm(0));
    Observations_ yobs(Jo.observations());
    std::unique_ptr<Departures_> ydep(Jo.newDepartures());
    yobs.obsvalues() -= ydep->depvalues(); 
    Log::test() << "H(xt): " << yobs << std::endl;
    
    // Perturb observations
    std::unique_ptr<Departures_> ypert(Jo.randomizeCovar());
    yobs += *ypert;
    Log::test() << "H(xt)+Eo: " << yobs << std::endl;

    // Write perturbed observations
    const std::vector<eckit::LocalConfiguration> obsConf =
      fullConfig.getSubConfigurations("cost_function.Jo");
    yobs.save(obsConf[0].getString("Observation.ObsData.obsvalue"));

    // Setup variables
    const std::vector<std::string> varNames = fullConfig.getStringVector("variables");
    JediVariables tmpVars(varNames);
    const Variables_ varsT(util::templatedVarsConf(tmpVars));

    // Create background perturbation
    Increment_ dx(resol, varsT, xx.state()[0].validTime());
    J->jb().jbState().covar().randomize(dx);

    // Write background perturbation
    const eckit::LocalConfiguration pertConf(fullConfig, "background perturbation");
    dx.write(pertConf);

    // Create configuration to write out first-guess departure
    const eckit::LocalConfiguration depConf(fullConfig, "first-guess departure");
    eckit::LocalConfiguration diagnostic;
    diagnostic.set("departures", depConf.getString("obsvalue"));
    eckit::LocalConfiguration evalConf;
    evalConf.set("diagnostics", diagnostic);

    // Re-evalute J aroung background instead of truth
    CtrlVariable_ xxTmp(xx);
    xxTmp.state()[0] += dx;
    J->evaluate(xxTmp, evalConf);

    return 0;
  }

  // -----------------------------------------------------------------------------
 private:
  std::string appname() const { return "oops::SyntheticData<" + MODEL::name() + ">"; }
  // -----------------------------------------------------------------------------
};

}  // namespace oops
