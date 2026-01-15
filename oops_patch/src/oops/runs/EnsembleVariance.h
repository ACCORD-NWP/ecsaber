/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Ensemble.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/Variables.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/ECUtilities.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class EnsembleVariance : public Application {
  using Ensemble_ = Ensemble<MODEL>;
  using Geometry_ = Geometry<MODEL>;
  using Increment_ = Increment<MODEL>;
  using Model_ = Model<MODEL>;
  using State_ = State<MODEL>;
  using Variables_ = Variables<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  EnsembleVariance() {}
  // -----------------------------------------------------------------------------
  virtual ~EnsembleVariance() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    //  Setup model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    //  Setup variables
    const eckit::LocalConfiguration varConfig(fullConfig, "variables");
    const Variables_ vars(varConfig);

    // Configurations for different ensembles
    std::vector<eckit::LocalConfiguration> confs =
        fullConfig.getSubConfigurations("ensembles");

    for (const auto& conf : confs) {
      // Set name
      Log::info() << "Set name" << std::endl;
      std::string name = conf.getString("name");

      // Read state
      Log::info() << "Read state" << std::endl;
      eckit::LocalConfiguration stateConf(conf, "state");
      const State_ state(resol, model, stateConf);

      // Read ensemble
      Log::info() << "Read ensemble" << std::endl;
      eckit::LocalConfiguration ensConfig(conf, "ensemble");
      const size_t nens = ensConfig.getInt("members");
      util::expandEnsembleTemplate(ensConfig, nens);
      util::setMPI(ensConfig, eckit::mpi::comm().size());
      Ensemble_ ensemble(state.validTime(), ensConfig);
      ensemble.linearize(state, resol);

      // Compute variance
      Increment_ variance = ensemble.variance();

      // Write variance
      Log::test() << "Write variance of " << name << ": " << variance.norm()
                        << std::endl;
      Log::test() << variance << std::endl;
      const eckit::LocalConfiguration varianceConf(conf, "variance");
      variance.write(varianceConf);
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsembleVariance<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops
