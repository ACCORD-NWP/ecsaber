/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

#include "saber/interpolation/GaussToCS.h"
#include "saber/interpolation/Rescaling.h"

using atlas::grid::detail::partitioner::TransPartitioner;
using atlas::grid::detail::partitioner::MatchingMeshPartitionerCubedSphere;

namespace saber {
namespace interpolation {

namespace {

int getHalo(const std::string & interpType) {
  // Assumes that only interpolation types with "cubic" in their name require a halo of 2.
  const std::string pattern = "cubic";
  if (interpType.find(pattern) != interpType.npos) {
    return 2;
  }
  return 1;
}

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid,
                             const int halo) {
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(halo));
}

// -----------------------------------------------------------------------------

auto createPointCloud(const atlas::Grid& grid,
                      const atlas::grid::Partitioner& partitioner) {
    const auto distribution = atlas::grid::Distribution(grid, partitioner);

    auto lonLats = std::vector<atlas::PointXY>{};
    auto idx = atlas::gidx_t{0};
    for (const auto& lonLat : grid.lonlat()) {
        if (distribution.partition(idx++) == atlas::mpi::rank()) {
            lonLats.emplace_back(lonLat.data());
        }
    }
    return std::make_unique<atlas::functionspace::PointCloud>(
                atlas::functionspace::PointCloud(lonLats));
}

// -----------------------------------------------------------------------------

auto createInverseInterpolation(const bool initializeInverseInterpolation,
                                const bool hasSinglePE,
                                const atlas::functionspace::NodeColumns & csFunctionSpace,
                                const atlas::Grid & gaussGrid,
                                const atlas::grid::Partitioner & gaussPartitioner) {
  CS2Gauss inverseInterpolation;

  if (!initializeInverseInterpolation || hasSinglePE) {
    return inverseInterpolation;
  }

  const auto config = atlas::util::Config("type", "cubedsphere");
  const atlas::grid::MatchingMeshPartitioner csMatchingPartitioner(csFunctionSpace.mesh(),
                                                                   config);
  inverseInterpolation.matchingPtcldFspace = createPointCloud(gaussGrid, csMatchingPartitioner);

  inverseInterpolation.targetPtcldFspace = createPointCloud(gaussGrid, gaussPartitioner);

  const auto interpConfig =  atlas::util::Config("type", "cubedsphere-bilinear") |
                             atlas::util::Config("halo_exchange", false);

  inverseInterpolation.interpolation = atlas::Interpolation(
              interpConfig, csFunctionSpace,
              *inverseInterpolation.matchingPtcldFspace);

  inverseInterpolation.redistribution = atlas::Redistribution(
              *inverseInterpolation.matchingPtcldFspace,
              *inverseInterpolation.targetPtcldFspace);

  return inverseInterpolation;
}

// -----------------------------------------------------------------------------

/* Direct interpolation from NodeColumn cubed-sphere FunctionSpace to
 * StructuredColumns is not possible on multiple PEs.
 * Here, the interpolation is done following this route:
 * CubedSphere FunctionSpace
 * -> matching PointCloud FunctionSpace using CS partitioner
 * -> target PointCloud FunctionSpace using Gauss grid partitioner
 * -> gauss StructuredColumns FunctionSpace.
 */
void inverseInterpolateMultiplePEs(
        const oops::JediVariables & variables,
        const CS2Gauss & inverseInterpolation,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  // Interpolate from source to matching PointCloud
  atlas::FieldSet matchingPtcldFset;
  for (const auto & var : variables) {
    auto matchingPtcldField = inverseInterpolation.matchingPtcldFspace->createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(srcFieldSet[var.name()].shape(1)));
    atlas::array::make_view<double, 2>(matchingPtcldField).assign(0.0);
    matchingPtcldFset.add(matchingPtcldField);
  }
  srcFieldSet.haloExchange();
  inverseInterpolation.interpolation.execute(srcFieldSet, matchingPtcldFset);

  // Redistribute from matching PointCloud to target PointCloud
  atlas::FieldSet targetPtcldFset;
  for (const auto & var : variables) {
    auto targetPtcldField = inverseInterpolation.targetPtcldFspace->createField<double>(
          atlas::option::name(var.name()) |
          atlas::option::levels(srcFieldSet[var.name()].shape(1)));
    targetPtcldFset.add(targetPtcldField);
  }
  inverseInterpolation.redistribution.execute(matchingPtcldFset, targetPtcldFset);

  // Copy from target PointCloud to gauss StructuredColumns
  for (const auto & var : variables) {
    atlas::Field gaussField = gaussFunctionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(srcFieldSet[var.name()].shape(1)));
    atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(targetPtcldFset[var.name()]));
    gaussField.set_dirty();  // atlas interpolation/redistribution above produces dirty halos
    newFieldSet.add(gaussField);
  }
}

// -----------------------------------------------------------------------------

void inverseInterpolateSinglePE(
        const oops::JediVariables & variables,
        const atlas::functionspace::NodeColumns & CSFunctionSpace,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const atlas::Grid & gaussGrid,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  const auto targetMesh = atlas::MeshGenerator("structured").generate(gaussGrid);
  const auto hybridFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);
  const auto scheme = atlas::util::Config("type", "cubedsphere-bilinear") |
                      atlas::util::Config("adjoint", "false");
  const auto interp = atlas::Interpolation(scheme, CSFunctionSpace, hybridFunctionSpace);

  atlas::FieldSet hybridFieldSet;
  for (const auto & var : variables) {
    atlas::Field hybridField = hybridFunctionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(srcFieldSet[var.name()].shape(1)));
    atlas::array::make_view<double, 2>(hybridField).assign(0.0);
    hybridFieldSet.add(hybridField);
  }

  srcFieldSet.haloExchange();
  interp.execute(srcFieldSet, hybridFieldSet);

  // Copy into StructuredColumns
  for (const auto & var : variables) {
    atlas::Field gaussField = gaussFunctionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(srcFieldSet[var.name()].shape(1)));
    atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(hybridFieldSet[var.name()]));
    gaussField.set_dirty();  // atlas interpolation above produces dirty halos
    newFieldSet.add(gaussField);
  }
}

// -----------------------------------------------------------------------------

Rescaling initRescaling(const GaussToCSParameters & params,
                        const eckit::mpi::Comm & comm,
                        const oops::JediVariables & activeVars,
                        const atlas::FunctionSpace & innerFspace,
                        const atlas::FunctionSpace & outerFspace,
                        const saber::interpolation::AtlasInterpWrapper & interp) {
  if (params.interpolationRescaling.value().is_initialized()) {
    const auto & conf = params.interpolationRescaling.value().value();
    if (conf.has("horizontal covariance profile file path")) {
      // Compute rescaling fields from input covariance profile
      return Rescaling(comm,
                       conf,
                       activeVars,
                       innerFspace,
                       outerFspace,
                       interp);
    } else if (conf.has("input atlas file")) {
      // Read rescaling field from atlas file
      const auto readConf = conf.getSubConfiguration("input atlas file");
      return Rescaling(comm,
                       readConf,
                       activeVars,
                       outerFspace);
    } else {
      throw eckit::UserError("Missing parameters to initialize a relevant Rescaling object",
                             Here());
    }
  } else {
    return Rescaling();
  }
}
}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussToCS> makerGaussToCS_("gauss to cubed-sphere-dual");

// -----------------------------------------------------------------------------
// Note that this is slower than this needs to be
// as we are need to create 2 grid objects (very slow)
// In the future it might make sense to include an atlas grid (if available) from
// the model in outerGeometryData.
GaussToCS::GaussToCS(const oops::GeometryData & outerGeometryData,
                     const oops::JediVariables & outerVars,
                     const eckit::Configuration & covarConf,
                     const Parameters_ & params,
                     const oops::FieldSet3D & xb,
                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerVars_(outerVars),
    activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
    CSFunctionSpace_(outerGeometryData.functionSpace()),
    gaussGrid_(params.gaussGridUid.value()),
    gaussFunctionSpace_(createGaussFunctionSpace(gaussGrid_,
                                                 getHalo(params.interpType.value()))),
    gaussPartitioner_(new TransPartitioner()),
    csgrid_(CSFunctionSpace_.mesh().grid()),
    interp_(gaussPartitioner_, gaussFunctionSpace_, csgrid_, CSFunctionSpace_,
            params.interpType.value()),
    inverseInterpolation_(createInverseInterpolation(
                              params.initializeInverseInterpolation.value(),
                              outerGeometryData.comm().size() == 1,
                              CSFunctionSpace_, gaussGrid_,
                              gaussPartitioner_)),
    rescaling_(initRescaling(params,
                             outerGeometryData.comm(),
                             activeVars_,
                             gaussFunctionSpace_,
                             CSFunctionSpace_,
                             interp_)),
    innerGeometryData_(gaussFunctionSpace_, outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(),
                       outerGeometryData.comm())

{
  oops::Log::trace() << classname() << "::GaussToCS starting" << std::endl;
  oops::Log::trace() << classname() << "::GaussToCS done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
    if (activeVars_.has(fieldname)) {
      gaussFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // On input: fieldset on gaussian mesh

  // Create fieldset on cubed-sphere mesh.

  atlas::FieldSet csFieldSet;
  for (const auto & var : activeVars_) {
    atlas::Field csField =
      CSFunctionSpace_.createField<double>(
          atlas::option::name(var.name()) |
          atlas::option::levels(gaussFieldSet[var.name()].shape(1)) |
          atlas::option::halo(1));
    atlas::array::make_view<double, 2>(csField).assign(0.0);
    csFieldSet.add(csField);
  }

  // Interpolate to cubed sphere
  interp_.execute(gaussFieldSet, csFieldSet);
  csFieldSet.set_dirty();  // atlas interpolation produces dirty halos

  for (const auto & var : activeVars_) {
    newFields.add(csFieldSet[var.name()]);
  }

  fieldSet.fieldSet() = newFields;

  // Apply optional rescaling
  rescaling_.execute(fieldSet);

  oops::Log::trace() << classname() << "::multiply done"
                     << fieldSet.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // Apply optional rescaling
  rescaling_.execute(fieldSet);

  // Create empty fieldSets
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet csFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
    if (activeVars_.has(fieldname)) {
      csFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // Create gauss fieldset
  atlas::FieldSet gaussFieldSet;
  for (const auto & var : activeVars_) {
    atlas::Field gaussField =
      gaussFunctionSpace_.createField<double>(atlas::option::name(var.name()) |
            atlas::option::levels(csFieldSet[var.name()].shape(1)) |
            atlas::option::halo(1));
    atlas::array::make_view<double, 2>(gaussField).assign(0.0);
    gaussFieldSet.add(gaussField);
  }

  // Adjoint of interpolation from gauss to dual cubed sphere
  interp_.executeAdjoint(gaussFieldSet, csFieldSet);

  for (const auto & var : activeVars_) {
    newFields.add(gaussFieldSet[var.name()]);
  }

  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  atlas::FieldSet newFieldSet = atlas::FieldSet();
  atlas::FieldSet srcFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       srcFieldSet.add(fieldSet[fieldname]);
     } else {
       newFieldSet.add(fieldSet[fieldname]);
     }
  }

  if (innerGeometryData_.comm().size() >= 2) {
    inverseInterpolateMultiplePEs(activeVars_, inverseInterpolation_,
                                  gaussFunctionSpace_,
                                  srcFieldSet, newFieldSet);
  } else {
    // A faster and more direct route is possible on a single PE
    inverseInterpolateSinglePE(activeVars_,
                               CSFunctionSpace_,
                               gaussFunctionSpace_,
                               gaussGrid_,
                               srcFieldSet, newFieldSet);
  }

  fieldSet.fieldSet() = newFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCS::generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                                  const oops::JediVariables & innerVars) const {
  oops::FieldSet3D fset(this->validTime(), innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCS::generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                                  const oops::JediVariables & outerVars) const {
  oops::FieldSet3D fset(this->validTime(), outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void GaussToCS::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
