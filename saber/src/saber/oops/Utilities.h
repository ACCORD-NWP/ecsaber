/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/range/adaptor/reversed.hpp>

#include "atlas/field.h"
#include "atlas/grid.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSets.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelData.h"
#include "oops/interface/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/ECUtilities.h"

namespace oops {
  class FieldSet3D;
}

namespace saber {

// -----------------------------------------------------------------------------

oops::JediVariables getActiveVars(const SaberBlockParametersBase & params,
                              const oops::JediVariables & defaultVars);

// -----------------------------------------------------------------------------

oops::JediVariables getUnionOfInnerActiveAndOuterVars(const SaberBlockParametersBase & params,
                                                  const oops::JediVariables & outerVars);

// -----------------------------------------------------------------------------

oops::JediVariables getInnerOnlyVars(const SaberBlockParametersBase & params,
                                 const oops::JediVariables & outerVars);

// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration &,
               const int &);

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

// -----------------------------------------------------------------------------

void expandEnsembleTemplate(eckit::LocalConfiguration &,
                            const size_t &);

// -----------------------------------------------------------------------------

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::JediVariables & vars);

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::JediVariables & varsToAllocate,
                           const oops::JediVariables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace);

// -----------------------------------------------------------------------------

template<typename MODEL>
oops::FieldSets readEnsemble(const oops::Geometry<MODEL> & geom,
                             const oops::JediVariables & modelvars,
                             const oops::State4D<MODEL> & xb,
                             const oops::State4D<MODEL> & fg,
                             const eckit::LocalConfiguration & inputConf,
                             const bool & iterativeEnsembleLoading,
                             eckit::LocalConfiguration & outputConf) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Fill output configuration and set ensemble size
  size_t nens = 0;
  size_t ensembleFound = 0;
  eckit::LocalConfiguration varConf;

  // Ensemble of states, perturbation using the mean
  std::vector<eckit::LocalConfiguration> ensembleConf;
  if (inputConf.has("ensemble")) {
    if (util::isVector(inputConf.getSubConfiguration("ensemble"))) {
      ensembleConf = inputConf.getSubConfigurations("ensemble");
    } else {
      ensembleConf.push_back(inputConf.getSubConfiguration("ensemble"));
    }
    nens = ensembleConf[0].getInt("members");
    for (auto & ensemble3DConf : ensembleConf) {
      expandEnsembleTemplate(ensemble3DConf, nens);
    }
    outputConf.set("ensemble", ensembleConf);
    varConf = ensembleConf[0];
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk
  std::vector<eckit::LocalConfiguration> ensemblePert;
  if (inputConf.has("ensemble pert")) {
    if (util::isVector(inputConf.getSubConfiguration("ensemble pert"))) {
      ensemblePert = inputConf.getSubConfigurations("ensemble pert");
    } else {
      ensemblePert.push_back(inputConf.getSubConfiguration("ensemble pert"));
    }
    nens = ensemblePert[0].getInt("members");
    for (auto & ensemble3DConf : ensembleConf) {
      expandEnsembleTemplate(ensemble3DConf, nens);
    }
    outputConf.set("ensemble", ensemblePert);
    varConf = ensemblePert[0];
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk on other geometry
  eckit::LocalConfiguration ensemblePertOtherGeom;
  if (inputConf.has("ensemble pert on other geometry")
          && inputConf.has("ensemble geometry")) {
    ensemblePertOtherGeom = inputConf.getSubConfiguration("ensemble pert on other geometry");

    // Bespoke validation, mimicking oops::IncrementEnsembleParameters<MODEL>
    ASSERT(ensemblePertOtherGeom.has("date"));
    ASSERT(ensemblePertOtherGeom.has("members from template")
           || ensemblePertOtherGeom.has("members"));
    ASSERT(!(ensemblePertOtherGeom.has("members from template")
           && ensemblePertOtherGeom.has("members")));

    if (ensemblePertOtherGeom.has("members")) {
      const auto members = ensemblePertOtherGeom.getSubConfigurations("members");
      nens = members.size();
      varConf = members[0];
    }

    if (ensemblePertOtherGeom.has("members from template")) {
      const auto members = ensemblePertOtherGeom.getSubConfiguration("members from template");
      ASSERT(members.has("nmembers"));
      ASSERT(members.has("pattern"));
      ASSERT(members.has("template"));
      nens = members.getInt("nmembers");
      varConf = members.getSubConfiguration("template");
    }

    outputConf.set("ensemble pert on other geometry", ensemblePertOtherGeom);
    outputConf.set("ensemble geometry",
                   inputConf.getSubConfiguration("ensemble geometry"));
    ++ensembleFound;
  }

  // Set ensemble size
  outputConf.set("ensemble size", nens);

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  oops::JediVariables vars(varConf.has("variables") ?
    oops::JediVariables{varConf.getStringVector("variables")} :
    modelvars);

  if (!iterativeEnsembleLoading) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble of states, perturbation using the mean
    if (ensembleConf.size() > 0) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;

      for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
        std::shared_ptr<oops::Ensemble<MODEL>> ens_k(new oops::Ensemble<MODEL>(xb[jsub].validTime(),
          ensembleConf[jsub]));
        ens_k->linearize(xb[jsub], geom);
        for (size_t ie = 0; ie < nens; ++ie) {
          (*ens_k)[ie] *= std::sqrt(static_cast<double>(nens-1));
        }
        oops::EnsemblesCollection<MODEL>::getInstance().put(xb[jsub].validTime(), ens_k);
      }
    }

    // Increment ensemble from increments on disk
    if (ensemblePert.size() > 0 && ensemblePertOtherGeom.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;

      for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
        std::shared_ptr<oops::Ensemble<MODEL>> ens_k(new oops::Ensemble<MODEL>(xb[jsub].validTime(),
          ensemblePert[jsub]));
        ens_k->build(xb[jsub], geom);
        ens_k->read();
        oops::EnsemblesCollection<MODEL>::getInstance().put(xb[jsub].validTime(), ens_k);
      }
    }

    if (ensembleConf.size() > 0 || ensemblePert.size() > 0) {
      // Transform Ensemble into FieldSets
      oops::Log::info() << "Info     : Transform Ensemble into FieldSets" << std::endl;
      std::vector<int> ensmems(nens);
      std::iota(ensmems.begin(), ensmems.end(), 0);
      oops::FieldSets fsetEns(xb, ensmems);
      return fsetEns;
    }

    // Increment ensemble from increments on disk on other geometry
    if (!ensemblePertOtherGeom.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments "
                        << "on disk on other geometry" << std::endl;
      const eckit::mpi::Comm & commGeom = eckit::mpi::comm();

      // Setup functionspace
      auto fspaceConf = inputConf.getSubConfiguration("ensemble geometry");
      atlas::Grid grid;
      atlas::grid::Partitioner partitioner;
      atlas::Mesh mesh;
      atlas::FunctionSpace fspace;
      atlas::FieldSet fieldset;
      util::setupFunctionSpace(commGeom, fspaceConf, grid, partitioner,
                               mesh, fspace, fieldset);

      // Setup variable sizes
      if (fspaceConf.has("groups")) {
        // Read level information from configuration
        const auto groups = fspaceConf.getSubConfigurations("groups");
        for (const auto & group : groups) {
          const int levels = group.getInt("levels");
          for (const auto & var : group.getStringVector("variables")) {
            if (vars.has(var)) {
              vars[var].setLevels(levels);
            }
          }
        }
        // Check all variables have been populated with level information
        for (const auto & var : vars) {
          if (var.getLevels() < 0) {
            std::stringstream ss;
            ss << "Invalid vertical level information for variable "
               << var << " in `ensemble geometry: groups`.";
            throw eckit::UserError(ss.str(), Here());
          }
        }
      } else {
        // Use level information from the model variables
        for (auto & var : vars) {
          var.setLevels(modelvars[var.name()].getLevels());
        }
      }

      // Read perturbations into oops::FieldSets
      oops::FieldSets fsetEns(fspace, vars, xb.times(),
                              ensemblePertOtherGeom,
                              commGeom, eckit::mpi::self());
      return fsetEns;
    }
  }
  // Return empty ensemble if none was returned before
  std::vector<util::DateTime> dates;
  std::vector<int> ensmems;
  oops::FieldSets fsetEns(dates, eckit::mpi::self(), ensmems, eckit::mpi::self());
  return fsetEns;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::JediVariables & vars,
                      const util::DateTime & date,
                      const eckit::LocalConfiguration & conf,
                      oops::FieldSet3D & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create variables
  oops::Variables<MODEL> varsT(templatedVarsConf(vars));

  // Create Increment
  oops::Increment<MODEL> dx(geom, varsT, date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset.shallowCopy(dx.increment().fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::JediVariables & vars,
                        const eckit::LocalConfiguration & conf,
                        const size_t & ie,
                        oops::FieldSet3D & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  if (conf.has("ensemble")) {
    // Ensemble of states passed as increments
    std::vector<eckit::LocalConfiguration> ensembleConf = 
      conf.getSubConfigurations("ensemble");
    std::vector<eckit::LocalConfiguration> membersConf =
      ensembleConf[0].getSubConfigurations("state");

    // Create variables
    oops::Variables<MODEL> varsT(templatedVarsConf(vars));

    // Read state as increment
    oops::Increment<MODEL> dx(geom, varsT, fset.validTime());
    dx.read(membersConf[ie]);

    // Copy FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from difference of two states
    std::vector<eckit::LocalConfiguration> ensembleConf
      = conf.getSubConfigurations("ensemble pert");
    std::vector<eckit::LocalConfiguration> membersConf =
      ensembleConf[0].getSubConfigurations("state");

    // Create variables
    oops::Variables<MODEL> varsT(templatedVarsConf(vars));

    // Read Increment
    oops::Increment<MODEL> dx(geom, varsT, fset.validTime());
    dx.read(membersConf[ie]);

    // Get FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert on other geometry")
          && conf.has("ensemble geometry")) {
    throw eckit::NotImplemented("readEnsembleMember not yet implemented for an"
                                "ensemble on a non-MODEL geometry", Here());
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound <= 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
