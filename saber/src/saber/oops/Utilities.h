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

#include "atlas/field.h"
#include "atlas/grid.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSets.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
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
#include "oops/util/ParallelFieldSetIO.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "oops/util/ECUtilities.h"

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

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi,
            const int & omp);

// -----------------------------------------------------------------------------

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::JediVariables & vars);

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::JediVariables & varsToAllocate,
                           const oops::JediVariables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace);

// -----------------------------------------------------------------------------

size_t getNensFromConfig(const eckit::Configuration & conf);

// -----------------------------------------------------------------------------

eckit::LocalConfiguration getEnsSubconfig(const eckit::Configuration & conf, size_t iens);

// -----------------------------------------------------------------------------

oops::FieldSets readEnsemble(const oops::GeometryData & geomData,
                             const oops::JediVariables & modelvars,
                             const std::vector<util::DateTime> & times,
                             const eckit::mpi::Comm & commTime,
                             const eckit::mpi::Comm & commEns,
                             const eckit::Configuration & inputConf);

// -----------------------------------------------------------------------------

template<typename MODEL>
oops::FieldSets readEnsemble(const oops::Geometry<MODEL> & geom,
                             const oops::JediVariables & modelvars,
                             const std::vector<util::DateTime> & times,
                             const eckit::mpi::Comm & commTime,
                             const eckit::mpi::Comm & commEns,
                             const eckit::Configuration & inputConf) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Fill output configuration and set ensemble size
  size_t ensembleFound = 0;
  eckit::LocalConfiguration varConf;

  // Ensemble of states, perturbation using the mean
  eckit::LocalConfiguration ensembleConf;
  if (inputConf.has("ensemble")) {
    ensembleConf = inputConf.getSubConfiguration("ensemble");
    varConf = getEnsSubconfig(ensembleConf, 0);
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk
  eckit::LocalConfiguration ensemblePertConf;
  if (inputConf.has("ensemble pert")) {
    ensemblePertConf = inputConf.getSubConfiguration("ensemble pert");
    varConf = getEnsSubconfig(ensemblePertConf, 0);
    ++ensembleFound;
  }

  // Increment ensemble from difference of two states
  eckit::LocalConfiguration ensembleBaseConf;
  eckit::LocalConfiguration ensemblePairsConf;
  if (inputConf.has("ensemble base") && inputConf.has("ensemble pairs")) {
    ensembleBaseConf = inputConf.getSubConfiguration("ensemble base");
    ensemblePairsConf = inputConf.getSubConfiguration("ensemble pairs");
    varConf = getEnsSubconfig(ensembleBaseConf, 0);
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk on other geometry
  eckit::LocalConfiguration ensemblePertOtherGeomConf;
  eckit::LocalConfiguration ensembleGeomConf;
  if (inputConf.has("ensemble pert on other geometry") && inputConf.has("ensemble geometry")) {
    ensemblePertOtherGeomConf = inputConf.getSubConfiguration("ensemble pert on other geometry");
    ensembleGeomConf = inputConf.getSubConfiguration("ensemble geometry");
    varConf = getEnsSubconfig(ensemblePertOtherGeomConf, 0);
    ++ensembleFound;
  }

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  oops::JediVariables vars(varConf.has("variables") ?
    oops::JediVariables{varConf.getStringVector("variables")} :
    modelvars);

  if (!inputConf.getBool("iterative ensemble loading", false)) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble of states, perturbation using the mean
    if (!ensembleConf.empty()) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;
      std::vector<eckit::LocalConfiguration> memConfs;
      std::vector<int> ensmems;
      for (size_t ie = 0; ie < getNensFromConfig(ensembleConf); ++ie) {
        memConfs.push_back(getEnsSubconfig(ensembleConf, ie));
        ensmems.push_back(ie);
      }
      oops::FieldSets fsetEns(geom, vars, times, ensmems, memConfs, true);
      return fsetEns;
    }

    // Increment ensemble from increments on disk
    if (!ensemblePertConf.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
      std::vector<eckit::LocalConfiguration> memConfs;
      std::vector<int> ensmems;
      for (size_t ie = 0; ie < getNensFromConfig(ensemblePertConf); ++ie) {
        memConfs.push_back(getEnsSubconfig(ensemblePertConf, ie));
        ensmems.push_back(ie);
      }
      oops::FieldSets fsetEns(geom, vars, times, ensmems, memConfs, false);
      return fsetEns;
    }

    // Increment ensemble from difference of two states
    if (!ensembleBaseConf.empty() && !ensemblePairsConf.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from difference of two states"
                        << std::endl;
      throw eckit::Exception("not implemented yet", Here());
    }

    // Increment ensemble from increments on disk on other geometry
    if (!ensemblePertOtherGeomConf.empty() && !ensembleGeomConf.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments "
                        << "on disk on other geometry" << std::endl;
      const eckit::mpi::Comm & commGeom = eckit::mpi::comm();

      // Setup functionspace
      atlas::Grid grid;
      atlas::grid::Partitioner partitioner;
      atlas::Mesh mesh;
      atlas::FunctionSpace fspace;
      atlas::FieldSet fieldset;
      util::setupFunctionSpace(commGeom, ensembleGeomConf, grid, partitioner,
                               mesh, fspace, fieldset);

      // Setup variable sizes
      if (ensembleGeomConf.has("groups")) {
        // Read level information from configuration
        const auto groups = ensembleGeomConf.getSubConfigurations("groups");
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

      if (varConf.has("parallel IO")) {
        util::ParallelFieldSetIO io(fspace,
                                    ensemblePertOtherGeomConf.getString("grid name"),
                                    util::ParallelFieldSetIO::Mode::Read);

        // Read perturbations into oops::FieldSets
        oops::FieldSets fsetEns(fspace, vars, io, times,
                                ensemblePertOtherGeomConf,
                                commGeom, eckit::mpi::self());

        return fsetEns;
      } else {
        oops::FieldSets fsetEns(fspace, vars, times,
                                ensemblePertOtherGeomConf,
                                commGeom, commTime);
        return fsetEns;
      }
    }
  }

  // Return empty ensemble if none was returned before
  std::vector<util::DateTime> dates;
  std::vector<int> ensmems;
  oops::FieldSets ensemble(dates, commTime, ensmems, commEns);
  return ensemble;
}

// -------------------------------------------------------------------------------------------------

/// Scale all ensemble members by 1/sqrt(denom) in place.
inline void scaleEnsemble(oops::FieldSets & ensemble, const double denom) {
  ensemble *= 1.0/std::sqrt(denom);
}

// -------------------------------------------------------------------------------------------------

/// Read ensemble and immediately scale members by 1/sqrt(scaling).
/// If "denominator for normalizing ensemble covariance" is present in inputConf, uses that value;
/// otherwise defaults to ensemble size - 1.
template<typename MODEL>
oops::FieldSets readAndScaleEnsemble(const oops::Geometry<MODEL> & geom,
                                     const oops::JediVariables & modelvars,
                                     const std::vector<util::DateTime> & times,
                                     const eckit::mpi::Comm & commTime,
                                     const eckit::mpi::Comm & commEns,
                                     const eckit::Configuration & inputConf) {
  oops::Log::trace() << "readAndScaleEnsemble starting" << std::endl;
  oops::FieldSets ensemble = readEnsemble(geom, modelvars, times, commTime, commEns, inputConf);
  double denom = static_cast<double>(ensemble.ens_size() - 1);
  if (inputConf.has("denominator for normalizing ensemble covariance")) {
    denom = inputConf.getDouble("denominator for normalizing ensemble covariance");
  }
  scaleEnsemble(ensemble, denom);
  oops::Log::trace() << "readAndScaleEnsemble done" << std::endl;
  return ensemble;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::JediVariables & vars,
                      const util::DateTime & date,
                      const eckit::Configuration & conf,
                      oops::FieldSet3D & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create Increment
  oops::Increment<MODEL> dx(geom, util::templatedVars<MODEL>(vars), date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset.shallowCopy(dx.increment().fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -----------------------------------------------------------------------------

void readEnsembleMember(const oops::GeometryData & geomData,
                        const oops::JediVariables & vars,
                        const eckit::Configuration & conf,
                        const size_t & ie,
                        oops::FieldSet3D & fset);

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::JediVariables & vars,
                        const eckit::Configuration & conf,
                        const size_t & ie,
                        oops::FieldSet3D & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  if (conf.has("ensemble")) {
    // Read state as increment
    eckit::LocalConfiguration memConf = getEnsSubconfig(
      conf.getSubConfiguration("ensemble"), ie);
    oops::Increment<MODEL> dx(geom, util::templatedVars<MODEL>(vars), fset.validTime());
    dx.read(memConf);

    // Copy FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from increments on disk
    eckit::LocalConfiguration memConf = getEnsSubconfig(
      conf.getSubConfiguration("ensemble pert"), ie);

    // Read Increment
    oops::Increment<MODEL> dx(geom, util::templatedVars<MODEL>(vars), fset.validTime());
    dx.read(memConf);

    // Get FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble base") && conf.has("ensemble pairs")) {
    throw eckit::NotImplemented("readEnsembleMember not yet implemented for a"
                                "difference of two states", Here());
  }

  if (conf.has("ensemble pert on other geometry")
          && conf.has("ensemble geometry")) {
    throw eckit::NotImplemented("readEnsembleMember not yet implemented for an"
                                "ensemble on a non-MODEL geometry", Here());
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound == 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

void cvToFset(const atlas::Field & cv,
              oops::FieldSet3D & fset,
              const size_t & offset,
              const oops::JediVariables & vars);

// -----------------------------------------------------------------------------

void fsetToCv(const oops::FieldSet3D & fset,
              atlas::Field & cv,
              const size_t & offset,
              const oops::JediVariables & vars);

// -----------------------------------------------------------------------------

}  // namespace saber
