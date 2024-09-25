/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------

double Fields::min(const Variables & vars) const {
  oops::Log::trace() << classname() << "::min starting" << std::endl;

  double zmin = std::numeric_limits<double>::max();
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars.variablesList()) {
    const atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      const auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zmin = std::min(zmin, view(jnode, jlevel));
          }
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(zmin, eckit::mpi::min());

  oops::Log::trace() << classname() << "::min done" << std::endl;
  return zmin;
}

// -----------------------------------------------------------------------------

double Fields::max(const Variables & vars) const {
  oops::Log::trace() << classname() << "::max starting" << std::endl;

  double zmax = -std::numeric_limits<double>::max();
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars.variablesList()) {
    const atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      const auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zmax = std::max(zmax, view(jnode, jlevel));
          }
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(zmax, eckit::mpi::max());

  oops::Log::trace() << classname() << "::max done" << std::endl;
  return zmax;
}

// -----------------------------------------------------------------------------

void Fields::interpolate(const Locations & locs,
                         GeoVaLs & gv) const {
  oops::Log::trace() << classname() << "::interpolate starting" << std::endl;

  if (locs.grid().size() > 0) {
    // Setup interpolation
    const auto & interpolation = setupObsInterpolation(locs);

    // Create observation fieldset
    atlas::FieldSet obsFieldSet;
    for (const auto & var : vars_.variablesList()) {
      atlas::Field obsField = interpolation->dstFspace().createField<double>(
        atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
      obsFieldSet.add(obsField);
    }

    // Copy FieldSet
    atlas::FieldSet fset = util::copyFieldSet(fset_);

    // Horizontal interpolation
    interpolation->execute(fset, obsFieldSet);

    // Vertical interpolation
    interpolation->executeVertical(obsFieldSet, gv.fieldSet());
  }

  oops::Log::trace() << classname() << "::interpolate done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::interpolateAD(const Locations & locs,
                           const GeoVaLs & gv) {
  oops::Log::trace() << classname() << "::interpolateAD starting" << std::endl;

  if (locs.grid().size() > 0) {
    // Setup interpolation
    const auto & interpolation = setupObsInterpolation(locs);

    // Create observation fieldset
    atlas::FieldSet obsFieldSet;
    for (const auto & var : vars_.variablesList()) {
      atlas::Field obsField = interpolation->dstFspace().createField<double>(
        atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
      obsFieldSet.add(obsField);
    }

    // Vertical interpolation
    interpolation->executeVerticalAdjoint(obsFieldSet, gv.fieldSet());

    // Horizontal interpolation
    interpolation->executeAdjoint(fset_, obsFieldSet);

    // Reduce duplicate points
    reduceDuplicatePoints();
  }

  oops::Log::trace() << classname() << "::interpolateAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::forceWith(const Fields & other,
                       const Variables & vars) {
  oops::Log::trace() << classname() << "::Fields forceWith" << std::endl;

  // Copy time
  time_ = other.time_;

  // Copy fields
  for (const auto & var : vars.variablesList()) {
    atlas::Field field = fset_[var];
    const atlas::Field fieldOther = other.fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) = viewOther(jnode, jlevel);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::forceWith done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::synchronizeFields() {
  oops::Log::trace() << classname() << "::synchronizeFields starting" << std::endl;

  // Set duplicate points to the same value
  resetDuplicatePoints();

  oops::Log::trace() << classname() << "::synchronizeFields done" << std::endl;
}

// -----------------------------------------------------------------------------

eckit::Stream & operator<<(eckit::Stream & s,
                           const Fields & rhs) {
  oops::Log::trace() << "operator<< starting" << std::endl;

  std::vector<double> vect;
  rhs.serialize(vect);
  for (auto & value : vect) {
    s << value;
  }

  oops::Log::trace() << "operator<< done" << std::endl;
  return s;
}

// -----------------------------------------------------------------------------

eckit::Stream & operator>>(eckit::Stream & s,
                           Fields & rhs) {
  oops::Log::trace() << "operator>> starting" << std::endl;

  std::vector<double> vect;
  vect.resize(rhs.serialSize());
  for (auto & value : vect) {
      s >> value;
  }
  size_t index = 0;
  rhs.deserialize(vect, index);

  oops::Log::trace() << "operator>> done" << std::endl;
  return s;
}

// -----------------------------------------------------------------------------

std::vector<Interpolation>::iterator Fields::setupObsInterpolation(const Locations & locs)
  const {
  oops::Log::trace() << classname() << "::setupObsInterpolation starting" << std::endl;

  // Get geometry UIDs (grid + "_" + paritioner)
  const std::string srcGeomUid = geom_->grid().uid() + "_" + geom_->partitioner().type();
  const std::string dstObsUid = locs.grid().uid() + "_" + geom_->partitioner().type();

  // Compare with existing UIDs
  for (auto it = interpolations().begin(); it != interpolations().end(); ++it) {
    if ((it->srcUid() == srcGeomUid) && (it->dstUid() == dstObsUid)) {
      oops::Log::trace() << classname() << "::setupGridInterpolation done" << std::endl;
      return it;
    }
  }

  // New grid
  const int nobsGlb = locs.grid().size();

  // Define partition
  std::vector<int> partition(nobsGlb);
  size_t joGlb = 0;
  for (size_t jt = 0; jt < geom_->getComm().size(); ++jt) {
    for (int joLoc = 0; joLoc < locs.size(jt); ++joLoc) {
      partition[joGlb] = jt;
      ++joGlb;
    }
  }

  // Create function space
  std::unique_ptr<atlas::FunctionSpace> fspace;

  if (nobsGlb > 3) {
    // Create observation distribution
    atlas::grid::Distribution distribution(geom_->getComm().size(), nobsGlb, &partition[0]);

    // Create observation mesh
    atlas::Mesh obsMesh = atlas::MeshGenerator("delaunay").generate(locs.grid(), distribution);

    // Create observation function space (NodeColumns)
    fspace.reset(new atlas::functionspace::NodeColumns(obsMesh));
  } else {
    // Create observation function space (PointCloud)
    fspace.reset(new atlas::functionspace::PointCloud(locs.grid()));
  }

  // Create horizontal interpolation
  Interpolation interpolation(geom_->interpolation(),
                              geom_->getComm(),
                              geom_->partitioner(),
                              geom_->functionSpace(),
                              srcGeomUid,
                              locs.grid(),
                              *fspace,
                              dstObsUid);

  // Interpolate vertical coordinate
  atlas::FieldSet fset;
  atlas::FieldSet fsetInterp;
  for (const auto & field : fset_) {
    const std::string vertCoordName = "vert_coord_"
      + std::to_string(geom_->groupIndex(field.name()));
    if (!fset.has(vertCoordName)) {
      fset.add(geom_->fields()[vertCoordName]);
      atlas::Field fieldInterp = fspace->createField<double>(
        atlas::option::name(vertCoordName) | atlas::option::levels(field.levels()));
      fsetInterp.add(fieldInterp);
    }
  }
  interpolation.execute(fset, fsetInterp);

  // Setup vertical interpolation
  for (const auto & var : vars_.variablesList()) {
    const std::string vertCoordName = "vert_coord_" + std::to_string(geom_->groupIndex(var));
    const auto vert_coordView = atlas::array::make_view<double, 2>(fsetInterp[vertCoordName]);
    std::vector<std::array<size_t, 2>> verStencil;
    std::vector<std::array<double, 2>> verWeights;
    std::vector<size_t> verStencilSize;
    verStencil.resize(locs.size());
    verWeights.resize(locs.size());
    verStencilSize.resize(locs.size());
    for (int jo = 0; jo < locs.size(); ++jo) {
      if (geom_->levels(var) == 1) {
        // No vertical interpolation
        verStencilSize[jo] = 0;
      } else {
        // Linear vertical interpolation
        const double z = locs[jo][2];
        double bottom = std::numeric_limits<double>::max();
        double top = -std::numeric_limits<double>::max();
        for (size_t k = 0; k < geom_->levels(var); ++k) {
          bottom = std::min(bottom, vert_coordView(jo, k));
          top = std::max(top, vert_coordView(jo, k));
        }
        ASSERT(z >= bottom);
        ASSERT(z <= top);
        double zinf = -std::numeric_limits<double>::max();
        double zsup = std::numeric_limits<double>::max();
        size_t kinf = 0;
        size_t ksup = std::numeric_limits<size_t>::max();
        for (size_t k = 0; k < geom_->levels(var); ++k) {
          const double level = vert_coordView(jo, k);
          if (level == z) {
            zinf = level;
            zsup = level;
            kinf = k;
            ksup = k;
          } else {
            if (z > level && zinf < level) {
              zinf = level;
              kinf = k;
            }
            if (z < level && zsup > level) {
              zsup = level;
              ksup = k;
            }
          }
        }
        if (kinf == ksup) {
          verStencil[jo][0] = kinf;
          verWeights[jo][0] = 1.0;
          verStencilSize[jo] = 1;
        } else {
          verStencil[jo][0] = kinf;
          verWeights[jo][0] = (zsup-z)/(zsup-zinf);
          verStencil[jo][1] = ksup;
          verWeights[jo][1] = (z-zinf)/(zsup-zinf);
          verStencilSize[jo] = 2;
        }
      }
    }
    interpolation.insertVerticalInterpolation(var,
                                              verStencil,
                                              verWeights,
                                              verStencilSize);
  }

  // Insert new interpolation
  interpolations().push_back(interpolation);

  oops::Log::trace() << classname() << "::setupObsInterpolation done" << std::endl;
  return std::prev(interpolations().end());
}

// -----------------------------------------------------------------------------

void Fields::reduceDuplicatePoints() {
  oops::Log::trace() << classname() << "::reduceDuplicatePoints starting" << std::endl;

  if (geom_->duplicatePoints()) {
    if (geom_->gridType() == "regular_lonlat") {
      // Deal with poles
      for (auto field_internal : fset_) {
        // Get local sums
        atlas::functionspace::StructuredColumns fs(field_internal.functionspace());
        atlas::StructuredGrid grid = fs.grid();
        auto view = atlas::array::make_view<double, 2>(field_internal);
        auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
        auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
        std::vector<double> north(field_internal.shape(1), 0.0);
        std::vector<double> south(field_internal.shape(1), 0.0);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (view_j(jnode) == 1) {
              for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                north[jlevel] += view(jnode, jlevel);
              }
            }
            if (view_j(jnode) == grid.ny()) {
              for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                south[jlevel] += view(jnode, jlevel);
              }
            }
          }
        }

        // Reduce
        geom_->getComm().allReduceInPlace(north.begin(), north.end(), eckit::mpi::sum());
        geom_->getComm().allReduceInPlace(south.begin(), south.end(), eckit::mpi::sum());

        // Copy value
        for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
          for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (view_i(jnode) == 1) {
              if (view_j(jnode) == 1) {
                for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                  view(jnode, jlevel) = north[jlevel];
                }
              }
              if (view_j(jnode) == grid.ny()) {
                for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                  view(jnode, jlevel) = south[jlevel];
                }
              }
            } else {
              if ((view_j(jnode) == 1) || (view_j(jnode) == grid.ny())) {
                for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                  view(jnode, jlevel) = 0.0;
                }
              }
            }
          }
        }
      }
    } else {
      throw eckit::NotImplemented("duplicate points not supported for this grid", Here());
    }
  }

  oops::Log::trace() << classname() << "::reduceDuplicatePoints done" << std::endl;
}

// -----------------------------------------------------------------------------
