/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "oops/base/FieldSet3D.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

oops::JediVariables getActiveVars(const SaberBlockParametersBase & params,
                              const oops::JediVariables & defaultVars) {
  oops::Log::trace() << "getActiveVars starting" << std::endl;
  oops::JediVariables activeVars_nomd;
  if (params.mandatoryActiveVars().size() == 0) {
    // No mandatory active variables for this block
    activeVars_nomd = params.activeVars.value().get_value_or(defaultVars);
  } else {
    // Block with mandatory active variables
    activeVars_nomd = params.activeVars.value().get_value_or(params.mandatoryActiveVars());
    ASSERT(params.mandatoryActiveVars() <= activeVars_nomd);
  }
  // Copy the variables that exist in defaultVars from defaultVars (they have metadata
  // associated with them)
  oops::JediVariables activeVars;
  for (auto & var : activeVars_nomd) {
    if (defaultVars.has(var.name())) {
      activeVars.push_back(defaultVars[var.name()]);
    } else {
      activeVars.push_back(var);
    }
  }
  return activeVars;
}

// -----------------------------------------------------------------------------

// Return inner variables as outer variables + inner active variables
// Can be used to help define innerVars_ member in SABER outer blocks
oops::JediVariables getUnionOfInnerActiveAndOuterVars(const SaberBlockParametersBase & params,
                                                  const oops::JediVariables & outerVars) {
  oops::JediVariables innerVars(outerVars);
  innerVars += params.activeInnerVars(outerVars);
  return innerVars;
}

// -----------------------------------------------------------------------------

// Return inner variables that are not outer variables
oops::JediVariables getInnerOnlyVars(const SaberBlockParametersBase & params,
                                 const oops::JediVariables & outerVars) {
  oops::JediVariables innerOnlyVars(getUnionOfInnerActiveAndOuterVars(params, outerVars));
  innerOnlyVars -= outerVars;
  return innerOnlyVars;
}

// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration & conf,
               const int & member) {
  oops::Log::trace() << "setMember starting" << std::endl;

  if (conf.has("member pattern")) {
    std::string memberPattern = conf.getString("member pattern");
    util::seekAndReplace(conf, memberPattern, std::to_string(member));
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

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::JediVariables & vars) {
  for (const auto& var : vars) {
    if (fset.has(var.name())) {
      throw eckit::UserError("Variable " + var.name() + " is already allocated in FieldSet.",
                             Here());
    }
  }
}

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::JediVariables & varsToAllocate,
                           const oops::JediVariables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace) {
  oops::Log::trace() << "allocateMissingFields starting" << std::endl;
  for (const auto& var : varsToAllocate) {
    if (!fset.has(var.name())) {
      oops::Log::info() << "Info     : Allocating " << var.name() << std::endl;
      auto field = functionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(varsWithLevels[var.name()].getLevels()));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      field.set_dirty(false);
      fset.add(field);
    }
  }
  oops::Log::trace() << "allocateMissingFields done" << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace saber
