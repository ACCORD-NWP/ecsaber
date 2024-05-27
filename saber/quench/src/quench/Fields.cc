/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Fields.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "quench/FieldsIOBase.h"
#include "quench/Geometry.h"
#include "quench/GeoVaLs.h"
#include "quench/Locations.h"

namespace quench {

// -----------------------------------------------------------------------------

static std::vector<quench::Interpolation> interpolations_;

// -----------------------------------------------------------------------------

std::vector<quench::Interpolation>& Fields::interpolations() {
  return interpolations_;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Geometry & geom,
               const Variables & vars,
               const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    fset_.add(field);
  }


  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  // Set fields to zero
  this->zero();

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other,
               const Geometry & geom)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Check number of levels
  for (const auto & var : vars_.variables()) {
    if (geom_->levels(var) != geom.levels(var)) {
      throw eckit::Exception("Different number of levels for variable " + var
        + ", cannot interpolate", Here());
    }
  }

  if (geom_->grid() == other.geom_->grid() && geom_->halo() == other.geom_->halo()) {
    // Copy fieldset
    fset_ = util::copyFieldSet(other.fset_);
  } else {
    // Setup interpolation
    const auto & interpolation = setupGridInterpolation(*other.geom_);

    // Create fieldset
    for (const auto & var : vars_.variables()) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
      fset_.add(field);
    }

    // Set interpolation type
    for (auto field : fset_) {
      field.metadata().set("interp_type", "default");
    }

    // Horizontal interpolation
    interpolation->execute(other.fset_, fset_);
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other,
               const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    fset_.add(field);
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (const auto & var : vars_.variables()) {
      atlas::Field field = fset_[var];
      const atlas::Field fieldOther = other.fset_[var];
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        const auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = viewOther(jnode, jlevel);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields and copy data
  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    const atlas::Field fieldOther = other.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewOther(jnode, jlevel);
        }
      }
    }
    fset_.add(field);
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::zero() {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
    }
  }
  fset_.set_dirty(false);

  oops::Log::trace() << "Fields::zero end" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::constantValue(const double & value) {
  oops::Log::trace() << classname() << "::constantValue starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) view(jnode, jlevel) = value;
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::constantValue end" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::constantValue(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::constantValue starting" << std::endl;

  for (const auto & group : config.getSubConfigurations("constant group-specific value")) {
    const std::vector<std::string> vars = group.getStringVector("variables");
    const double value = group.getDouble("constant value");
    for (const auto & var : vars_.variables()) {
      if (std::find(vars.begin(), vars.end(), var) != vars.end()) {
        atlas::Field field = fset_[var];
        const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
        const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
        if (field.rank() == 2) {
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(0.0);
          for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              if (gmaskView(jnode, jlevel) == 1) view(jnode, jlevel) = value;
            }
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::constantValue end" << std::endl;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewRhs(jnode, jlevel);
        }
      }
      field.set_dirty(fieldRhs.dirty());
    }
  }
  time_ = rhs.time_;

  oops::Log::trace() << classname() << "::operator= end" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator+=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  // Right-hand side fieldset
  atlas::FieldSet fsetRhs;
  if (geom_->grid() == rhs.geom_->grid() && geom_->halo() == rhs.geom_->halo()) {
    // Same geometry
    fsetRhs = util::shareFields(rhs.fset_);
  } else {
    // Interpolate
    const Fields rhsInterp(rhs, *geom_);

    // Copy fieldset
    fsetRhs = util::copyFieldSet(rhsInterp.fset_);
  }

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldRhs = fsetRhs[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) += viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator-=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator-= starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) -= viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator*=(const double & zz) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) *= zz;
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void Fields::axpy(const double & zz,
                  const Fields & rhs) {
  oops::Log::trace() << classname() << "::axpy starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) += zz * viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

double Fields::dot_product_with(const Fields & fld2) const {
  oops::Log::trace() << classname() << "::dot_product_with starting" << std::endl;

  double zz = 0;
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    atlas::Field field1 = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field field2 = fld2.fset_[var];
    if (field1.rank() == 2) {
      auto view1 = atlas::array::make_view<double, 2>(field1);
      auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zz += view1(jnode, jlevel)*view2(jnode, jlevel);
          }
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  oops::Log::trace() << classname() << "::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void Fields::schur_product_with(const Fields & dx) {
  oops::Log::trace() << classname() << "::schur_product_with starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldDx = dx.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewDx = atlas::array::make_view<double, 2>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) *= viewDx(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldDx.dirty());
    }
  }

  oops::Log::trace() << classname() << "::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::random() {
  oops::Log::trace() << classname() << "::random starting" << std::endl;

  fset_.clear();
  for (size_t groupIndex = 0; groupIndex < geom_->groups(); ++groupIndex) {
    // Mask and ghost points fields
    const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());

    // Total size
    size_t n = 0;
    std::vector<std::string> groupVars;
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        groupVars.push_back(var);
      }
    }
    for (const auto & var : groupVars) {
      atlas::Field field = fset_[var];
      if (field.rank() == 2) {
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) ++n;
          }
        }
      }
    }
    geom_->getComm().allReduceInPlace(n, eckit::mpi::sum());

    // Local masks
    atlas::FieldSet localMasks;
    localMasks.add(geom_->fields()[gmaskName]);
    localMasks.add(geom_->functionSpace().ghost());

    // Global masks
    atlas::FieldSet globalMasks;
    atlas::Field gmaskGlobal = geom_->functionSpace().createField<int>(
      atlas::option::name(gmaskName) | atlas::option::levels(geom_->levels(groupIndex))
      | atlas::option::global());
    globalMasks.add(gmaskGlobal);
    atlas::Field ghostGlobal = geom_->functionSpace().createField<int>(atlas::option::name("ghost")
     | atlas::option::global());
    globalMasks.add(ghostGlobal);

    // Global data
    atlas::FieldSet globalData;
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels(var)) | atlas::option::global());
        globalData.add(field);
      }
    }

    // Gather masks on main processor
    if (geom_->functionSpace().type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      fs.gather(localMasks, globalMasks);
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      // NodeColumns
      if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
        fs.gather(localMasks, globalMasks);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(geom_->functionSpace());
        fs.gather(localMasks, globalMasks);
      }
    } else {
      throw eckit::NotImplemented(geom_->functionSpace().type() +
        " function space not supported yet", Here());
    }

    if (geom_->getComm().rank() == 0) {
      // Random vector
      util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);

      // Copy random values
      n = 0;
      const auto ghostView = atlas::array::make_view<int, 1>(globalMasks["ghost"]);
      for (const auto & var : vars_.variables()) {
        if (geom_->groupIndex(var) == groupIndex) {
          atlas::Field field = globalData[var];
          const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
          const auto gmaskView = atlas::array::make_view<int, 2>(globalMasks[gmaskName]);
          if (field.rank() == 2) {
            auto view = atlas::array::make_view<double, 2>(field);
            for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
                  view(jnode, jlevel) = rand_vec[n];
                  ++n;
                }
              }
            }
          }
        }
      }
    }

    // Local data
    atlas::FieldSet localData;
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels(var)));
        localData.add(field);
      }
    }

    // Scatter data from main processor
    if (geom_->functionSpace().type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      fs.scatter(globalData, localData);
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      // NodeColumns
      if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
        fs.scatter(globalData, localData);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(geom_->functionSpace());
        fs.scatter(globalData, localData);
      }
    } else {
      throw eckit::NotImplemented(geom_->functionSpace().type() +
        " function space not supported yet", Here());
    }

    // Copy data
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        fset_.add(localData[var]);
      }
    }
  }

  fset_.set_dirty();  // code is too complicated, mark dirty to be safe

  oops::Log::trace() << "Fields::random done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::dirac(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::dirac starting" << std::endl;

  // Get dirac specifications
  std::vector<double> lon = config.getDoubleVector("lon");
  std::vector<double> lat = config.getDoubleVector("lat");
  std::vector<atlas::idx_t> level = config.getIntVector("level");
  std::vector<std::string> variable = config.getStringVector("variable");

  // Check sizes
  if (lon.size() != lat.size()) throw eckit::UserError("Inconsistent dirac specification size",
    Here());
  if (lon.size() != level.size()) throw eckit::UserError("Inconsistent dirac specification size",
    Here());
  if (lon.size() != variable.size()) throw eckit::UserError("Inconsistent dirac specification size",
    Here());

  // Build KDTree for each MPI task
  atlas::util::IndexKDTree search;
  search.reserve(geom_->functionSpace().size());
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  const auto lonlatView = atlas::array::make_view<double, 2>(geom_->functionSpace().lonlat());
  atlas::idx_t n = 0;
  for (atlas::idx_t jnode = 0; jnode < geom_->functionSpace().size(); ++jnode) {
    if (ghostView(jnode) == 0) {
      atlas::PointLonLat pointLonLat(lonlatView(jnode, 0), lonlatView(jnode, 1));
      pointLonLat.normalise();
      atlas::PointXY point(pointLonLat);
      search.insert(point, n++);
    }
  }
  search.build();

  // Set fields to zero
  this->zero();

  // Set dirac points
  for (size_t jdir = 0; jdir < lon.size(); ++jdir) {
    // Get field
    atlas::Field field = fset_[variable[jdir]];

    // Find MPI task
    atlas::PointLonLat pointLonLat(lon[jdir], lat[jdir]);
    pointLonLat.normalise();

    // Search nearest neighbor
    size_t index = std::numeric_limits<size_t>::max();
    double distance = std::numeric_limits<double>::max();
    if (geom_->functionSpace().size() > 0) {
      atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);
      index = neighbor[0].payload();
      distance = neighbor[0].distance();
    }
    std::vector<double> distances(geom_->getComm().size());
    geom_->getComm().gather(distance, distances, 0);

    // Find local task
    size_t localTask(-1);
    if (geom_->getComm().rank() == 0) {
      localTask = std::distance(std::begin(distances), std::min_element(std::begin(distances),
        std::end(distances)));
    }
    geom_->getComm().broadcast(localTask, 0);

    if (geom_->getComm().rank() == localTask) {
      // Add Dirac impulse
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        view(index, level[jdir]-1) = 1.0;
      }
    }
  }

  oops::Log::trace() << classname() << "::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::diff(const Fields & x1,
                  const Fields & x2) {
  oops::Log::trace() << classname() << "::diff starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldx1 = x1.fset_[var];
    atlas::Field fieldx2 = x2.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewx1 = atlas::array::make_view<double, 2>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 2>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) = viewx1(jnode, jlevel)-viewx2(jnode, jlevel);
          }
        }
      }
      field.set_dirty(fieldx1.dirty() || fieldx2.dirty());
    }
  }

  oops::Log::trace() << classname() << "::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::toFieldSet starting" << std::endl;

  // Share internal fieldset
  fset.clear();
  fset = util::shareFields(fset_);
  for (auto field_external : fset) {
    field_external.metadata().set("interp_type", "default");
  }
  fset.set_dirty(false);

  oops::Log::trace() << classname() << "::toFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::trace() << classname() << "::fromFieldSet starting" << std::endl;

  // Check input fieldset
  ASSERT(!fset.empty());

  // Reset internal fieldset
  fset_.clear();
  fset_ = util::shareFields(fset);

  // Reset variables
  vars_ = Variables(fset_.field_names());

  if (geom_->gridType() == "regular_lonlat") {
    // Reset poles points
    for (auto field_internal : fset_) {
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
          if ((view_j(jnode) == 1)  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              north[jlevel] = view(jnode, jlevel);
            }
          }
          if ((view_j(jnode) == grid.ny())  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              south[jlevel] = view(jnode, jlevel);
            }
          }
        }
      }
      geom_->getComm().allReduceInPlace(north.begin(), north.end(), eckit::mpi::sum());
      geom_->getComm().allReduceInPlace(south.begin(), south.end(), eckit::mpi::sum());
      for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
          atlas::idx_t jnode = fs.index(i, j);
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
        }
      }
    }
  }
  oops::Log::trace() << classname() << "::fromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::read(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Get IO format
  const std::string ioFormat = config.getString("format", "default");

  // Set FieldsIO
  std::unique_ptr<FieldsIOBase> fieldsIO(FieldsIOFactory::create(ioFormat));

  // Update variables names
  std::vector<std::string> variableNames;
  for (const auto & var : vars_.variables()) {
    std::string newVar = var;
    for (const auto & item : geom_->alias()) {
      if (item.getString("in code") == var) {
        newVar = item.getString("in file");
      }
    }
    variableNames.push_back(newVar);
  }

  // Create variables
  eckit::LocalConfiguration varConf;
  varConf.set("variables", variableNames);
  Variables vars(varConf);

  // Read fields
  fieldsIO->read(*geom_, vars, config, fset_);

  // Rename fields
  for (auto & field : fset_) {
    for (const auto & item : geom_->alias()) {
      if (item.getString("in file") == field.name()) {
        field.rename(item.getString("in code"));
      }
    }
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  fset_.set_dirty();

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::write(const eckit::Configuration & config) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Get IO format
  const std::string ioFormat = config.getString("format", "default");

  // Set FieldsIO
  std::unique_ptr<FieldsIOBase> fieldsIO(FieldsIOFactory::create(ioFormat));

  // Copy fieldset
  atlas::FieldSet fset = util::copyFieldSet(fset_);

  // Rename fields
  for (auto & field : fset) {
    for (const auto & item : geom_->alias()) {
      if (item.getString("in code") == field.name()) {
        field.rename(item.getString("in file"));
      }
    }
  }

  // Write fields
  fieldsIO->write(*geom_, config, fset);

  if (geom_->mesh().generated() && config.getBool("write gmsh", false)) {
    // GMSH file path
    std::string gmshfilepath = config.getString("filepath");;
    gmshfilepath.append(".msh");
    oops::Log::info() << "Info     : Writing file: " << gmshfilepath << std::endl;

    // GMSH configuration
    const auto gmshConfig =
    atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", true) |
    atlas::util::Config("info", true);
    atlas::output::Gmsh gmsh(gmshfilepath, gmshConfig);

     // Write GMSH
    gmsh.write(geom_->mesh());
    gmsh.write(fset, fset[0].functionspace());
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

double Fields::norm() const {
  oops::Log::trace() << classname() << "::norm" << std::endl;
  return util::normFieldSet(fset_, vars_.variables(), geom_->getComm());
}

// -----------------------------------------------------------------------------

void Fields::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os.setf(std::ios::scientific);
  os.precision(std::numeric_limits<double>::digits10+1);
  os << std::endl;
  os << *geom_;
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix << "Fields:";
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    os << std::endl;
    double zz = 0.0;
    atlas::Field field = fset_[var];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zz += view(jnode, jlevel)*view(jnode, jlevel);
          }
        }
      }
    }
    geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
    zz = sqrt(zz);
    os << prefix << "  " << var << ": " << zz;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t Fields::serialSize() const {
  oops::Log::trace() << classname() << "::serialSize starting" << std::endl;

  size_t nn = 0;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      nn += field.shape(0)*field.shape(1);
    }
  }

  oops::Log::trace() << classname() << "::serialSize done" << std::endl;
  return nn;
}

// -----------------------------------------------------------------------------

void Fields::serialize(std::vector<double> & vect)  const {
  oops::Log::trace() << classname() << "::serialize starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    const atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      const auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          vect.push_back(view(jnode, jlevel));
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::deserialize(const std::vector<double> & vect,
                         size_t & index) {
  oops::Log::trace() << classname() << "::deserialize starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = vect[index];
          ++index;
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<Interpolation>::iterator Fields::setupGridInterpolation(const Geometry & srcGeom)
  const {
  oops::Log::trace() << classname() << "::setupGridInterpolation starting" << std::endl;

  // Compare with exisiting UIDs
  for (auto it = interpolations().begin(); it != interpolations().end(); ++it) {
    if ((it->srcUid() == srcGeom.grid().uid()) && (it->dstUid() == geom_->grid().uid())) {
      oops::Log::trace() << classname() << "::setupGridInterpolation done" << std::endl;
      return it;
    }
  }

  // Create interpolation
  Interpolation interpolation(geom_->getComm(), srcGeom.partitioner(), srcGeom.functionSpace(),
    geom_->grid(), geom_->functionSpace());

  // Insert new interpolation
  interpolations().push_back(interpolation);

  oops::Log::trace() << classname() << "::setupGridInterpolation done" << std::endl;
  return std::prev(interpolations().end());
}

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

    // Horizontal interpolation
    interpolation->execute(fset_, obsFieldSet);

    // Vertical interpolation
    for (const auto & var : vars_.variablesList()) {
      const auto obsView = atlas::array::make_view<double, 2>(obsFieldSet[var]);
      auto gvView = atlas::array::make_view<double, 1>(gv.fieldSet()[var]);
      gvView.assign(0.0);
      const std::vector<interp::InterpElement> & verInterp = interpolation->verticalInterpolation(var);
      for (int jo = 0; jo < locs.size(); ++jo) {
        for (const auto & operation : verInterp[jo].operations()) {
          gvView(jo) += operation.second*obsView(jo, operation.first);
        }
      }
    }
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

    // Vertical interpolation adjoint
    for (const auto & var : vars_.variablesList()) {
      const auto gvView = atlas::array::make_view<double, 1>(gv.fieldSet()[var]);
      auto obsView = atlas::array::make_view<double, 2>(obsFieldSet[var]);
      obsView.assign(0.0);
      const std::vector<interp::InterpElement> & verInterp = interpolation->verticalInterpolation(var);
      for (int jo = 0; jo < locs.size(); ++jo) {
        for (const auto & operation : verInterp[jo].operations()) {
          obsView(jo, operation.first) += operation.second*gvView(jo);
        }
      }
    }

    // Horizontal interpolation
    interpolation->executeAdjoint(fset_, obsFieldSet);
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

  // Check internal fieldset
  ASSERT(!fset_.empty());

  if (geom_->gridType() == "regular_lonlat") {
    // Reset poles points
    for (auto field_internal : fset_) {
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
          if ((view_j(jnode) == 1)  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              north[jlevel] = view(jnode, jlevel);
            }
          }
          if ((view_j(jnode) == grid.ny())  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              south[jlevel] = view(jnode, jlevel);
            }
          }
        }
      }
      geom_->getComm().allReduceInPlace(north.begin(), north.end(), eckit::mpi::sum());
      geom_->getComm().allReduceInPlace(south.begin(), south.end(), eckit::mpi::sum());
      for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
          atlas::idx_t jnode = fs.index(i, j);
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
        }
      }
    }
  }

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

  // Compare with exisiting UIDs
  for (auto it = interpolations().begin(); it != interpolations().end(); ++it) {
    if ((it->srcUid() == geom_->grid().uid()) && (it->dstUid() == locs.grid().uid())) {
      oops::Log::trace() << classname() << "::setupObsInterpolation done" << std::endl;
      return it;
    }
  }

  // New grid
  const int nobsGlb = locs.grid().size();

  // Define partition
  std::vector<int> partition(nobsGlb);
  size_t joGlb = 0;
  for (size_t jt = 0; jt < geom_->getComm().size(); ++jt) {
    for (size_t joLoc = 0; joLoc < locs.size(jt); ++joLoc) {
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
  Interpolation interpolation(geom_->getComm(), geom_->partitioner(), geom_->functionSpace(),
    locs.grid(), *fspace);

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
    std::vector<interp::InterpElement> verInterp;
    for (int jo = 0; jo < locs.size(); ++jo) {
      std::vector<std::pair<size_t, double>> operations;
      if (geom_->levels(var) == 1) {
        // No vertical interpolation
        operations.push_back(std::make_pair(0, 1.0));
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
          operations.push_back(std::make_pair(kinf, 1.0));
        } else {
          operations.push_back(std::make_pair(kinf, (zsup-z)/(zsup-zinf)));
          operations.push_back(std::make_pair(ksup, (z-zinf)/(zsup-zinf)));
        }
      }
      verInterp.push_back(interp::InterpElement(operations));
    }
    interpolation.insertVerticalInterpolation(var, verInterp);
  }

  // Insert new interpolation
  interpolations().push_back(interpolation);

  oops::Log::trace() << classname() << "::setupObsInterpolation done" << std::endl;
  return std::prev(interpolations().end());
}

// -----------------------------------------------------------------------------

}  // namespace quench
