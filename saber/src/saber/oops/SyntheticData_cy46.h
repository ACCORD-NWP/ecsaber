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

#include "oops/base/Departures.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/PostProcessor.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsErrorCovariance.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/ECUtilities.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SyntheticData : public oops::Application {
  using CovarianceBase_ = oops::ModelSpaceCovarianceBase<MODEL>;
  using CovarianceFactory_ = oops::CovarianceFactory<MODEL>;
  using Departures_ = oops::Departures<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using ModelAux_ = oops::ModelAuxControl<MODEL>;
  using ObsAuxCtrl_ = oops::ObsAuxControl<MODEL>;
  using ObsErrorCovariance_ = oops::ObsErrorCovariance<MODEL>;
  using Observations_ = oops::Observations<MODEL>;
  using ObsSpace_ = oops::ObservationSpace<MODEL>;
  using ObsOperator_ = oops::ObsOperator<MODEL>;
  using ObsVector_ = oops::ObsVector<MODEL>;
  using State_ = oops::State<MODEL>;
  using Variables_ = oops::Variables<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  SyntheticData() {
    oops::instantiateObsErrorFactory<MODEL>();
    saber::instantiateCovarFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~SyntheticData() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    //  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "assimilation_window");
    const util::Duration winlen(windowConf.getString("window_length"));
    const util::DateTime winbgn(windowConf.getString("window_begin"));
    const util::DateTime winend(winbgn + winlen);
    oops::Log::info() << "Observation window is:" << windowConf << std::endl;

    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    //  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    //  Setup truth state
    const eckit::LocalConfiguration truthConfig(fullConfig, "truth");
    oops::Log::info() << "Truth configuration is:" << truthConfig << std::endl;
    State_ xt(resol, model, truthConfig);
    oops::Log::test() << "Truth state: " << xt.norm() << std::endl;

    //  Setup augmented state
    const ModelAux_ moderr(resol, truthConfig);

    //  Setup forecast outputs
    oops::PostProcessor<State_> postt;
    oops::PostProcessor<State_> postb;
    
    //  Setup observations bias
    const eckit::LocalConfiguration biasConf(fullConfig, "ObsBias");
    const ObsAuxCtrl_ ybias(biasConf);

    //  Setup observations
    std::vector<boost::shared_ptr<oops::Observer<MODEL, State_> > > pobst;
    std::vector<boost::shared_ptr<oops::Observer<MODEL, State_> > > pobsb;
    std::vector<std::unique_ptr<ObsSpace_>> obsdb;
    const std::vector<eckit::LocalConfiguration> obsConfs = fullConfig.getSubConfigurations("Observations");
    const size_t nobs = obsConfs.size();

    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      // Observations configuration
      eckit::LocalConfiguration obsConf(obsConfs[jobs], "Observation");
      oops::Log::debug() << "Observation configuration is:" << obsConf << std::endl;
      obsdb.emplace_back(new ObsSpace_(obsConf, winbgn, winend));
      const ObsOperator_ hop(*obsdb[jobs], obsConf);
      const Observations_ yy(*obsdb[jobs]);
      boost::shared_ptr<oops::Observer<MODEL, State_> >
        ppt(new oops::Observer<MODEL, State_>(*obsdb[jobs], hop, yy, ybias));
      boost::shared_ptr<oops::Observer<MODEL, State_> >
        ppb(new oops::Observer<MODEL, State_>(*obsdb[jobs], hop, yy, ybias));
      postt.enrollProcessor(ppt);
      postb.enrollProcessor(ppb);
      pobst.push_back(ppt);
      pobsb.push_back(ppb);
    }

    //  Setup variables
    const std::vector<std::string> varNames = fullConfig.getStringVector("variables");
    oops::JediVariables tmpVars(varNames);
    const Variables_ varsT(templatedVarsConf(tmpVars));

    //  Setup background error covariance
    const eckit::LocalConfiguration bMatConf(fullConfig, "Covariance");
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(bMatConf, resol, varsT, xt));
    Bmat->linearize(xt, resol, bMatConf);

    //  Create background perturbation
    Increment_ dx(resol, varsT, xt.validTime());
    Bmat->randomize(dx);
    State_ xb(xt);
    xb += dx;

    //  Write background perturbation
    const eckit::LocalConfiguration pertConf(fullConfig, "background perturbation");
    dx.write(pertConf);

    //  Compute H(xt)
    model.forecast(xt, moderr, winlen, postt);
    oops::Log::info() << "HofX: Finished observation computation." << std::endl;
    oops::Log::test() << "Final state (xt): " << xt.norm() << std::endl;

    //  Compute H(xb)
    model.forecast(xb, moderr, winlen, postb);
    oops::Log::info() << "HofX: Finished observation computation." << std::endl;
    oops::Log::test() << "Final state (xb): " << xb.norm() << std::endl;

    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      // Get observations from the truth
      boost::scoped_ptr<Observations_> yobst(pobst[jobs]->release());
      oops::Log::test() << "H(xt): " << *yobst << std::endl;
    
      // Setup R matrix
      eckit::LocalConfiguration rMatConf(obsConfs[jobs], "Covariance");
      ObsErrorCovariance_ Rmat(*obsdb[jobs], rMatConf);
      Rmat.linearize(*yobst);
      
      // Perturb observations
      Departures_ ydep(*obsdb[jobs]);
      Rmat.randomize(ydep);
      *yobst += ydep;
      oops::Log::test() << "H(xt)+Eo: " << *yobst << std::endl;

      // Write perturbed observations
      yobst->save("obsvalue");

      // Get observations from the background
      boost::scoped_ptr<Observations_> yobsb(pobsb[jobs]->release());
      oops::Log::test() << "H(xb): " << *yobsb << std::endl;

      // Compute (H(xt)+Eo)-H(xb)
      yobst->obsvalues() -= yobsb->obsvalues();
      yobst->save("fg_depar@body");
      oops::Log::test() << "(H(xt)+Eo)-H(xb): " << *yobst << std::endl;
    }

    return 0;
  }

  // -----------------------------------------------------------------------------
 private:
  std::string appname() const { return "saber::SyntheticData<" + MODEL::name() + ">"; }
  // -----------------------------------------------------------------------------
};

}  // namespace saber
