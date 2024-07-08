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

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Departures.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObservationSpaces.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class MakeObs : public oops::Application {
  using Departures_ = oops::Departures<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using ModelAux_ = oops::ModelAuxControl<MODEL>;
  using ObsAuxCtrls_ = oops::ObsAuxControls<MODEL>;
  using Observations_ = oops::Observations<MODEL>;
  using ObsSpaces_ = oops::ObservationSpaces<MODEL>;
  using ObsOperators_ = oops::ObsOperators<MODEL>;
  using State_ = oops::State<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  MakeObs() { oops::instantiateObsErrorFactory<MODEL>(); }
  // -----------------------------------------------------------------------------
  virtual ~MakeObs() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    //  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig,
                                               "Assimilation Window");
    const util::DateTime bgn(windowConf.getString("Begin"));
    const util::DateTime end(windowConf.getString("End"));
    const util::Duration fclen(end - bgn);
    oops::Log::info() << "Observation window is:" << windowConf << std::endl;

    //  Setup resolution
    const eckit::LocalConfiguration geomConfig(fullConfig, "Geometry");
    const Geometry_ geom(geomConfig);

    //  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(geom, modelConfig);

    //  Setup initial "true" state
    const eckit::LocalConfiguration initialConfig(fullConfig,
                                                  "Initial Condition");
    oops::Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    State_ xx(geom, model, initialConfig);
    oops::Log::test() << "Initial state: " << xx.norm() << std::endl;

    //  Setup augmented state
    ModelAux_ moderr(geom, model, initialConfig);

    //  Setup forecast outputs
    oops::PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("prints", prtConf);
    post.enrollProcessor(new oops::StateInfo<State_>("fc", prtConf));

    //  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    oops::Log::info() << "Observation configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, geom, bgn, end);
    ObsOperators_ hop(obspace);

    //  Setup observation bias
    ObsAuxCtrls_ ybias(obspace, obsconf);

    std::shared_ptr<oops::Observer<MODEL, State_> > pobs(
        new oops::Observer<MODEL, State_>(obspace, hop, ybias));
    post.enrollProcessor(pobs);

    //  Run forecast and release observations
    model.forecast(xx, moderr, fclen, post);
    oops::Log::test() << "Final state: " << xx.norm() << std::endl;
    oops::Log::info() << "MakeObs: Finished observation generation." << std::endl;
    std::unique_ptr<Observations_> yobs(pobs->release());

    // Setup empty departures
    Departures_ ydep(obspace);
    std::vector<Departures_> ypertVec;

    // Setup observation errors
    oops::ObsErrors<MODEL> matR(obspace);

    // Generate perturbations, first step
    matR.randomize(ydep);
    *yobs += ydep;

    // Member replacement
    std::string pattern = obsconf.getString("pattern");
    size_t zpad = obsconf.getInt("zero padding", 0);

    // Save perturbed observations
    for (std::size_t jj = 0; jj < yobs->size(); ++jj) {
      oops::Log::test() << "Generated observation: " << (*yobs)[jj] << std::endl;
    }
    eckit::LocalConfiguration obsconfControl(obsconf);
    util::seekAndReplace(obsconfControl, pattern, 0, zpad);
    yobs->save(obsconfControl);

    // Generate perturbations, second step
    size_t members = obsconf.getInt("members", 0);

    // Observation perturbation inflation factor
    double obspert = obsconf.getDouble("obspert", 1.0);

    if (members > 1) {
      // Generate and store observation perturbations
      for (size_t ie = 0; ie < members; ++ie) {
        matR.randomize(ydep);
        ydep *= obspert;
        ypertVec.push_back(ydep);
      }

      // Compute and subtract observation perturbations mean
      Departures_ ypertMean(obspace);
      for (size_t ie = 0; ie < members; ++ie) {
        ypertMean += ypertVec[ie];
      }
      ypertMean *= 1.0 / static_cast<double>(members);
      for (size_t ie = 0; ie < members; ++ie) {
        ypertVec[ie] -= ypertMean;
      }
    }

    for (size_t ie = 0; ie < members; ++ie) {
      // Add perturbation
      *yobs += ypertVec[ie];

      //  Save observations
      for (std::size_t jj = 0; jj < yobs->size(); ++jj) {
        oops::Log::test() << "Generated perturbed observation: " << (*yobs)[jj] << std::endl;
      }
      eckit::LocalConfiguration obsconfMember(obsconf);
      util::seekAndReplace(obsconfMember, pattern, ie+1, zpad);
      yobs->save(obsconfMember);

      // Remove perturbation
      *yobs -= ypertVec[ie];
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const { return "oops::MakeObs<" + MODEL::name() + ">"; }
  // -----------------------------------------------------------------------------
};

}  // namespace saber
