/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Geometry.h"

#include <netcdf.h>

#include <cmath>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "src/Variables.h"

#include "oops/generic/gc99.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "util/Logger.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace quench {

// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config)
  : comm_(eckit::mpi::comm()), groups_() {
  oops::Log::trace() << "Geometry::Geometry starting" << std::endl;

  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Setup atlas geometric data structures
  atlas::FieldSet fieldsetOwnedMask;
  util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh_, functionSpace_,
    fieldsetOwnedMask);

  halo_ = config.getInt("halo", 0);
  const eckit::LocalConfiguration gridParams(config, "grid");
  gridType_ = gridParams.getString("type", "no_type");
  regionalGrid_ = (gridType_ == "regional");  // grid.type() does not report if grid is regional

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : config.getSubConfigurations("groups")) {
    // Use this group index for all the group variables
    for (const auto & var : groupParams.getStringVector("variables")) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        throw eckit::UserError("Same variable present in distinct groups", Here());
      } else {
        groupIndex_[var] = groupIndex;
      }
    }

    // Define group
    groupData group;

    // Number of levels
    group.levels_ = groupParams.getInt("levels", 1);

    // Corresponding level for 2D variables (first or last)
    group.lev2d_ = groupParams.getString("lev2d", "first");

    // Vertical coordinate
    if (groupParams.has("vert_coord")) {
      std::vector<double> vert_coordParams = groupParams.getDoubleVector("vert_coord");
      if (vert_coordParams.size() != group.levels_) {
        throw eckit::UserError("Wrong number of levels in the user-specified vertical coordinate",
          Here());
      }
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        group.vert_coord_.push_back(vert_coordParams[jlevel]);
      }
    } else {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        group.vert_coord_.push_back(static_cast<double>(jlevel+1));
      }
    }

    // Default mask, set to 1 (true)
    atlas::Field gmask = functionSpace_.createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(group.levels_));
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    maskView.assign(1);

    // Specific mask
    std::string maskType = groupParams.getString("mask type", "none");
    if (maskType == "none") {
      // No mask
    } else if (maskType == "sea") {
      // Read sea mask
      readSeaMask(groupParams.getString("mask path", "../quench/data/landsea.nc"),
        group.levels_, group.lev2d_, gmask);
    } else {
      throw eckit::UserError("Wrong mask type", Here());
    }

    // Fill geometry fields
    group.fields_ = atlas::FieldSet();

    // Add owned points mask -- this mask does not depend on the group so was precomputed
    group.fields_->add(fieldsetOwnedMask.field("owned"));

    if (regionalGrid_) {
      // 2D indices
      const atlas::functionspace::StructuredColumns fs(functionSpace_);
      group.fields_->add(fs.index_i());
      group.fields_->add(fs.index_j());

      // Area
      atlas::Field area = functionSpace_.createField<double>(
        atlas::option::name("area") | atlas::option::levels(1));
      auto areaView = atlas::array::make_view<double, 2>(area);
      const atlas::StructuredGrid grid = fs.grid();
      const auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      const auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (atlas::idx_t jnode = 0; jnode < area.shape(0); ++jnode) {
        // Initialization
        int i = view_i(jnode)-1;
        int j = view_j(jnode)-1;
        double dist_i = 0.0;
        double dist_j = 0.0;

        // i-direction component
        if (i == 0) {
          dist_i = atlas::util::Earth().distance(grid.lonlat(i, j), grid.lonlat(i+1, j));
        } else if (i == grid.nx(j)-1) {
          dist_i = atlas::util::Earth().distance(grid.lonlat(i-1, j), grid.lonlat(i, j));
        } else {
          dist_i = 0.5*atlas::util::Earth().distance(grid.lonlat(i-1, j), grid.lonlat(i+1, j));
        }

        // j-direction component
        if (j == 0) {
          dist_j = atlas::util::Earth().distance(grid.lonlat(i, j), grid.lonlat(i, j+1));
        } else if (j == grid.ny()-1) {
          dist_j = atlas::util::Earth().distance(grid.lonlat(i, j-1), grid.lonlat(i, j));
        } else {
          dist_j = 0.5*atlas::util::Earth().distance(grid.lonlat(i, j-1), grid.lonlat(i, j+1));
       }

        // Local scale
        areaView(jnode, 0) = dist_i*dist_j;
      }
      group.fields_->add(area);
    }

    // Vertical coordinate
    atlas::Field vert_coord = functionSpace_.createField<double>(
      atlas::option::name("vert_coord") | atlas::option::levels(group.levels_));
    auto vert_coordView = atlas::array::make_view<double, 2>(vert_coord);
    for (atlas::idx_t jnode = 0; jnode < vert_coord.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        vert_coordView(jnode, jlevel) = group.vert_coord_[jlevel];
      }
    }

    // Add orography (mountain) on bottom level
    if (groupParams.has("orography")) {
      const eckit::LocalConfiguration orographyParams(groupParams, "orography");
      const atlas::PointLonLat topPoint({orographyParams.getDouble("top longitude"),
        orographyParams.getDouble("top latitude")});
      const double delta = (group.levels_ == 1) ? 1.0 :
        group.vert_coord_[group.levels_-2]-group.vert_coord_[group.levels_-1];
      const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
      for (atlas::idx_t jnode = 0; jnode < lonlatView.shape(0); ++jnode) {
        const atlas::PointLonLat xPoint({lonlatView(jnode, 0),
          orographyParams.getDouble("top latitude")});
        const atlas::PointLonLat yPoint({orographyParams.getDouble("top longitude"),
          lonlatView(jnode, 1)});
        double dxNorm = atlas::util::Earth().distance(xPoint, topPoint)
          /orographyParams.getDouble("zonal length");
        double dyNorm = atlas::util::Earth().distance(yPoint, topPoint)
          /orographyParams.getDouble("meridional length");
        double distNorm = std::sqrt(dxNorm*dxNorm+dyNorm*dyNorm);
        double orography = delta*orographyParams.getDouble("height")*oops::gc99(distNorm);
        vert_coordView(jnode, group.vert_coord_[group.levels_-1]) += orography;
      }
    }
    group.fields_->add(vert_coord);

    // Geographical mask
    group.fields_->add(gmask);

    // Mask size
    group.gmaskSize_ = 0.0;
    size_t domainSize = 0.0;
    auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
    for (atlas::idx_t jnode = 0; jnode < gmask.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < gmask.shape(1); ++jlevel) {
        if (ghostView(jnode) == 0) {
          if (maskView(jnode, jlevel) == 1) {
            group.gmaskSize_ += 1.0;
          }
          domainSize++;
        }
      }
    }
    comm_.allReduceInPlace(group.gmaskSize_, eckit::mpi::sum());
    comm_.allReduceInPlace(domainSize, eckit::mpi::sum());
    if (domainSize > 0) {
      group.gmaskSize_ = group.gmaskSize_/static_cast<double>(domainSize);
    }

    // Save group
    groups_.push_back(group);

    // Increment group index
    groupIndex++;
  }

  // Print summary
  this->print(oops::Log::info());

  oops::Log::trace() << "Geometry::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other)
  : comm_(other.comm_), halo_(other.halo_), grid_(other.grid_), gridType_(other.gridType_),
  regionalGrid_(other.regionalGrid_), partitioner_(other.partitioner_), mesh_(other.mesh_),
  groupIndex_(other.groupIndex_)  {
  oops::Log::trace() << "Geometry::Geometry starting" << std::endl;

  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
    } else {
      // Other NodeColumns
      functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    }
  } else if (other.functionSpace_.type() == "PointCloud") {
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported",
      Here());
  } else {
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported yet",
      Here());
  }

  // Copy groups
  for (size_t groupIndex = 0; groupIndex < other.groups_.size(); ++groupIndex) {
    // Define group
    groupData group;

    // Copy number of levels
    group.levels_ = other.groups_[groupIndex].levels_;

    // Copy corresponding level for 2D variables (first or last)
    group.lev2d_ = other.groups_[groupIndex].lev2d_;

    // Copy vertical coordinate
    group.vert_coord_ = other.groups_[groupIndex].vert_coord_;

    // Copy geometry fields
    group.fields_ = atlas::FieldSet();
    if (other.groups_[groupIndex].fields_.has("owned")) {
      group.fields_->add(other.groups_[groupIndex].fields_["owned"]);
    }
    group.fields_->add(other.groups_[groupIndex].fields_["vert_coord"]);
    group.fields_->add(other.groups_[groupIndex].fields_["gmask"]);

    // Copy mask size
    group.gmaskSize_ = other.groups_[groupIndex].gmaskSize_;

    // Save group
    groups_.push_back(group);
  }

  oops::Log::trace() << "Geometry::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t Geometry::levels(const std::string & var) const {
  oops::Log::trace() << "Geometry::levels starting" << std::endl;

  if (groupIndex_.count(var) == 0) {
    throw eckit::Exception("Variable " + var + " not found in groupIndex_", Here());
  }

  oops::Log::trace() << "Geometry::levels done" << std::endl;
  return groups_[groupIndex_.at(var)].levels_;
}

// -----------------------------------------------------------------------------

size_t Geometry::groupIndex(const std::string & var) const {
  oops::Log::trace() << "Geometry::groupIndex starting" << std::endl;

  if (groupIndex_.count(var) == 0) {
    throw eckit::Exception("Variable " + var + " not found in groupIndex_", Here());
  }

  oops::Log::trace() << "Geometry::groupIndex done" << std::endl;
  return groupIndex_.at(var);
}

// -----------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const Variables & vars) const {
  oops::Log::trace() << "Geometry::variableSizes starting" << std::endl;

  std::vector<size_t> sizes;
  for (const auto & var : vars.variablesList()) {
    sizes.push_back(levels(var));
  }

  oops::Log::trace() << "Geometry::variableSizes done" << std::endl;
  return sizes;
}

// -----------------------------------------------------------------------------

void Geometry::latlon(std::vector<double> & lats,
                      std::vector<double> & lons,
                      const bool includeHaloForRealLife) const {
  oops::Log::trace() << "Geometry::latlon starting" << std::endl;

  const auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  const auto ghost = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // TODO(Algo): Remove/fix the hack below when GeometryData local KD tree needs
  // to be set up correctly (e.g. when UnstructuredInterpolator is used).
  // For now never include halo in the latlon output because halo points from
  // some atlas grids (e.g. gaussian) can have unrealistic latitudes (e.g. more
  // than 90 degrees) and those latitudes can't be handled by KD trees.
  // Global KD trees created in GeometryData are used for communication and
  // don't need halo information.
  // Local KD trees in GeometryData need halo information but aren't used unless
  // UnstructuredInterpolator is used.
  bool includeHalo = false;
  const size_t npts = functionSpace_.size();
  const size_t nptsReturned = [&]() {
    if (includeHalo && comm_.size() > 1) {
      return npts;
    } else {
      size_t result = 0;
      for (atlas::idx_t i = 0; i < ghost.shape(0); ++i) {
        if (ghost(i) == 0) {
          result++;
        }
      }
      return result;
    }
  }();

  lats.resize(nptsReturned);
  lons.resize(nptsReturned);

  size_t count = 0;
  for (size_t jj = 0; jj < npts; ++jj) {
    // copy owned points, i.e. points with ghost==?
    if (ghost(jj) == 0 || (includeHalo && comm_.size() > 1)) {
      lats[count] = lonlat(jj, 1);
      lons[count] = lonlat(jj, 0);
      if (lons[count] < 0.0) lons[count] += 360.0;
      count++;
    }
  }
  ASSERT(count == nptsReturned);

  oops::Log::trace() << "Geometry::latlon done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  oops::Log::trace() << "Geometry::print starting" << std::endl;

  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix <<  "Quench geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (regionalGrid_) {
    os << prefix << "Regional grid detected" << std::endl;
  }
  if (partitioner_) {
    os << prefix << "Partitioner:" << std::endl;
    os << prefix << "- type: " << partitioner_.type() << std::endl;
  }
  os << prefix << "Function space:" << std::endl;
  os << prefix << "- type: " << functionSpace_.type() << std::endl;
  os << prefix << "- halo: " << halo_ << std::endl;
  os << prefix << "Groups: " << std::endl;
  for (size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
    os << prefix << "- Group " << groupIndex << ":" << std::endl;
    os << prefix << "  Vertical levels: " << std::endl;
    os << prefix << "  - number: " << levels(groupIndex) << std::endl;
    os << prefix << "  - vert_coord: " << groups_[groupIndex].vert_coord_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }

  oops::Log::trace() << "Geometry::print done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::readSeaMask(const std::string & maskPath,
                           const size_t & levels,
                           const std::string & lev2d,
                           atlas::Field & gmask) const {
  oops::Log::trace() << "Geometry::readSeaMask starting" << std::endl;

  // Lon/lat sizes
  size_t nlon = 0;
  size_t nlat = 0;

  // NetCDF IDs
  int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(maskPath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

    // Get lon/lat sizes
    if ((retval = nc_inq_dimid(ncid, "lon", &nlon_id))) ERR(retval);
    if ((retval = nc_inq_dimid(ncid, "lat", &nlat_id))) ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nlon_id, &nlon))) ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nlat_id, &nlat))) ERR(retval);
  }

  // Broadcast lon/lat sizes
  comm_.broadcast(nlon, 0);
  comm_.broadcast(nlat, 0);

  // Coordinates and land-sea mask
  std::vector<double> lon(nlon);
  std::vector<double> lat(nlat);
  std::vector<int> lsm(nlat*nlon);

  if (comm_.rank() == 0) {
    // Get lon/lat
    if ((retval = nc_inq_varid(ncid, "lon", &lon_id))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, "lat", &lat_id))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, "LSMASK", &lsm_id))) ERR(retval);

    // Read data
    float zlon[nlon][1];
    float zlat[nlat][1];
    uint8_t zlsm[nlat][nlon];
    if ((retval = nc_get_var_float(ncid, lon_id, &zlon[0][0]))) ERR(retval);
    if ((retval = nc_get_var_float(ncid, lat_id, &zlat[0][0]))) ERR(retval);
    if ((retval = nc_get_var_ubyte(ncid, lsm_id, &zlsm[0][0]))) ERR(retval);

    // Copy data
    for (size_t ilon = 0; ilon < nlon; ++ilon) {
      lon[ilon] = zlon[ilon][0];
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
      lat[ilat] = zlat[ilat][0];
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
     for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lsm[ilat*nlon+ilon] = static_cast<int>(zlsm[ilat][ilon]);
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);
  }

  // Broadcast coordinates and land-sea mask
  comm_.broadcast(lon.begin(), lon.end(), 0);
  comm_.broadcast(lat.begin(), lat.end(), 0);
  comm_.broadcast(lsm.begin(), lsm.end(), 0);

  // Build KD-tree
  atlas::Geometry geometry(atlas::util::Earth::radius());
  atlas::util::IndexKDTree2D search(geometry);
  search.reserve(nlat*nlon);
  std::vector<double> lon2d;
  std::vector<double> lat2d;
  std::vector<size_t> payload2d;
  int jnode = 0;
  for (size_t ilat = 0; ilat < nlat; ++ilat) {
    for (size_t ilon = 0; ilon < nlon; ++ilon) {
      lon2d.push_back(lon[ilon]);
      lat2d.push_back(lat[ilat]);
      payload2d.push_back(jnode);
      ++jnode;
    }
  }
  search.build(lon2d, lat2d, payload2d);

  // Ghost points
  atlas::Field ghost = functionSpace_.ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);

  if (functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(functionSpace_);
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Find nearest neighbor
        size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
          lonlatView(jnode, 1)}).payload();

        // Ocean points for all levels
        for (size_t jlevel = 0; jlevel < levels; ++jlevel) {
          if (lsm[nn] == 0) {
             maskView(jnode, jlevel) = 1;
           } else {
             maskView(jnode, jlevel) = 0;
           }
         }

        // Ocean + small islands for:
        // - the first level of 3D fields,
        // - the 2D fields if lev2d = "first"
        if (lsm[nn] == 3) {
          if ((levels > 1) || (lev2d == "first")) {
            maskView(jnode, 0) = 1;
          }
        }
      }
    }
  } else {
    throw eckit::NotImplemented("Sea mask not supported for " + functionSpace_.type() + " yet",
      Here());
  }

  oops::Log::trace() << "Geometry::readSeaMask done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
