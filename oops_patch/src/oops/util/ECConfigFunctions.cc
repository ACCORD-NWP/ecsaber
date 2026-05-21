/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "oops/util/ECConfigFunctions.h"

#include <algorithm>
#include <string>
#include <vector>

#include "oops/util/ConfigFunctions.h"
#include "util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration & conf,
               const size_t & member) {
  oops::Log::trace() << "setMember starting" << std::endl;

  if (conf.has("member pattern")) {
    std::string memberPattern = conf.getString("member pattern");
    seekAndReplace(conf, memberPattern, std::to_string(member));
  } else {
    conf.set("member", member);
  }

  oops::Log::trace() << "setMember done" << std::endl;
}

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi) {
  oops::Log::trace() << "setMPI starting" << std::endl;

  if (conf.has("mpi pattern")) {
    std::string mpiPattern = conf.getString("mpi pattern");
    util::seekAndReplace(conf, mpiPattern, std::to_string(mpi));
  }

  oops::Log::trace() << "setMPI done" << std::endl;
}

// -----------------------------------------------------------------------------

void expandEnsembleTemplate(eckit::LocalConfiguration & conf,
                            const size_t & nens) {
  oops::Log::trace() << "expandEnsembleTemplate starting" << std::endl;

  if (conf.has("state from template")) {
    eckit::LocalConfiguration templateConf(conf, "state from template");
    std::vector<eckit::LocalConfiguration> stateConf;
    for (size_t ie = 0; ie < nens; ++ie) {
      // Get correct index
      size_t count = templateConf.getInt("start", 1);
      std::vector<int> except = templateConf.getIntVector("except", {});
      for (size_t jj = 0; jj <= ie; ++jj) {
        // Check for excluded members
        while (std::count(except.begin(), except.end(), count)) {
          count += 1;
        }

        // Update counter
        if (jj < ie) count += 1;
      }

      // Replace pattern recursively in the configuration
      eckit::LocalConfiguration memberConf(templateConf, "template");
      std::string pattern = templateConf.getString("pattern");
      size_t zpad = templateConf.getInt("zero padding", 0);
      util::seekAndReplace(memberConf, pattern, count, zpad);

      // Add member
      stateConf.push_back(memberConf);
    }

    // Add 3D ensemble
    conf.set("state", stateConf);
  }

  oops::Log::trace() << "expandEnsembleTemplate done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util
