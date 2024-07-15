/*
 * (C) Copyright 2011 ECMWF
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/runs/Application.h"
#include "oops/runs/Variational.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "util/MPIWrapper.h"

namespace saber {

template <typename MODEL>
class AssimEnsemble : public oops::Application {
 public:
  // -----------------------------------------------------------------------------
  AssimEnsemble() {}

  // -----------------------------------------------------------------------------
  virtual ~AssimEnsemble() {}

  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration& fullConfig) const {
    std::vector<eckit::LocalConfiguration> assimConfig;
    size_t nb_assim;
    if (fullConfig.has("assim member")) {
      // Explicit configurations
      fullConfig.get("assim member", assimConfig);
      nb_assim = assimConfig.size();
    } else {
      // Configurations template
      eckit::LocalConfiguration templateConfig(fullConfig, "assim member template");

      // Number of members
      nb_assim = fullConfig.getInt("assim size");

      // Zero padding
      const size_t zpad = fullConfig.getInt("assim zero padding");

      // Pattern
      const std::string pattern = fullConfig.getString("assim pattern");

      for (size_t ie = 0; ie < nb_assim; ++ie) {
        // Update template
        eckit::LocalConfiguration config(templateConfig);
        util::seekAndReplace(config, pattern, ie+1, zpad);
        assimConfig.push_back(config);
      }
    }

    if (fullConfig.getBool("assim sequential", false)) {
      // Loop over cases
      for (size_t iassim = 0; iassim < nb_assim; ++iassim) {
        oops::Log::info() << "Performing minim for member:" << iassim << " out of "
                          << nb_assim << std::endl;
        oops::Log::info() << "Configuration: " << assimConfig[iassim] << std::endl;

        oops::Variational<MODEL> var;
        var.execute(assimConfig[iassim]);
      }
      return 0;
    } else {
      // Run depending on MPI wrapper rank
      util::MPIWrapper& mpiw = util::MPIWrapper::Instance();
      size_t iassim = mpiw.my_rank();
      ASSERT(iassim < nb_assim);
      oops::Log::info() << "Performing minim for member:" << iassim << " out of "
                        << nb_assim << std::endl;
      oops::Log::info() << "Configuration: " << assimConfig[iassim] << std::endl;

      oops::Variational<MODEL> var;
      return var.execute(assimConfig[iassim]);
    }
  }

  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "saber::AssimEnsemble<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace saber
