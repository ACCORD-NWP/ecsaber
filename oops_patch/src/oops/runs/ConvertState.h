/*
 * (C) Copyright 2018-2021 UCAR
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <omp.h>

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/Geometry.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

template <typename MODEL> class ConvertStateStatesParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ConvertStateStatesParameters, Parameters)
  typedef State<MODEL> State_;

 public:
  RequiredParameter<eckit::LocalConfiguration> input{"input", this};
  RequiredParameter<eckit::LocalConfiguration> output{"output", this};
};

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

/// Options taken by the ConvertState application.
template <typename MODEL> class ConvertStateParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ConvertStateParameters, Parameters)

 public:
  /// Input Geometry parameters.
  RequiredParameter<eckit::LocalConfiguration> inputGeometry{"input geometry", this};

  /// Output Geometry parameters.
  RequiredParameter<eckit::LocalConfiguration> outputGeometry{"output geometry", this};

  /// Model
  RequiredParameter<eckit::LocalConfiguration> model{"model", this};

  /// Variable change parameters (and option to do inverse).
  OptionalParameter<eckit::LocalConfiguration> varChange{"variable change", this};

  /// States to be converted
  RequiredParameter<std::vector<ConvertStateStatesParameters<MODEL>>> states{"states", this};
};

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL>               Geometry_;
  typedef Model<MODEL>                  Model_;
  typedef State<MODEL>                  State_;
  typedef VariableChange<MODEL>         VariableChange_;
  typedef ConvertStateParameters<MODEL> ConvertStateParameters_;
  typedef ConvertStateStatesParameters<MODEL> ConvertStateStatesParameters_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertState() : Application() {}
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    ConvertStateParameters_ params;
    params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup resolution for input and output
    const Geometry_ resol1(params.inputGeometry);
    const Geometry_ resol2(params.outputGeometry);

    // Setup model
    const Model_ model(resol1, params.model);

    // Setup change of variable
    std::unique_ptr<VariableChange_> vc;
    oops::JediVariables varout;
    bool inverse = false;
    if (params.varChange.value() != boost::none) {
      eckit::LocalConfiguration chconf(params.varChange.value().value());
      if (chconf.has("output variables")) {
        vc.reset(new VariableChange_(chconf, resol2));
        varout = JediVariables(chconf, "output variables");
        inverse = chconf.getBool("do inverse", false);
      }
    }

//  List of input and output states
    const int nstates = params.states.value().size();

//  Loop over states
    for (int jm = 0; jm < nstates; ++jm) {
//    Read current state parameters
      const ConvertStateStatesParameters_ stateParams = params.states.value()[jm];

//    Print output
      Log::info() << "Converting state " << jm+1 << " of " << nstates << std::endl;

//    Read state
      State_ xxi(resol1, model, stateParams.input.value());
      Log::test() << "Input state: " << xxi << std::endl;

//    Copy and change resolution
      State_ xx(resol2, xxi);

//    Variable transform(s)
      if (vc) {
          // Create variable change
        oops::JediVariables varin = xx.state().variables();
        if (inverse) {
          vc->changeVarInverse(xx, varout);
        } else {
          vc->changeVar(xx, varout);
        }
        Log::test() << "Variable transform: " << *vc << std::endl;
        Log::test() << "Variable change from: " << varin << std::endl;
        Log::test() << "Variable change to: " << varout << std::endl;
        Log::test() << "State after variable transform: " << xx << std::endl;
      }

//    Write state
      eckit::LocalConfiguration outconf(stateParams.toConfiguration(), "output");
      xx.write(outconf);

      Log::test() << "Output state: " << xx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
