/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/GeometryData.h"

#include <algorithm>

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/Point.h"

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace detail {
using atlas::Point3;
bool pointOutsideCircumcircle(const Point3 & point, const Point3 & a,
                              const Point3 & b, const Point3 & c) {
  atlas::PointXYZ center = Point3::cross(a, b) + Point3::cross(b, c) + Point3::cross(c, a);
  const double normA = Point3::norm(a);
  const double normCenter = Point3::norm(center);

  // Protect against the degenerate case where points coincide and the triangle has zero area.
  // This can occur for grids with singularities. If the triangle has zero area, then the
  // circumcircle is not well defined anyway, so we arbitrarily choose to return false. By
  // returning early, we avoid a div-by-zero when re-scaling `center` below.
  // Note: 0.5*norm(center) gives the triangle area; take sqrt to go back to linear units.
  if (sqrt(0.5 * normCenter) < 1e-14 * normA) {
    return false;
  }

  // Scale center to lie on sphere
  center *= normA / normCenter;

  // Check if point is outside circumcircle by comparing chord distances from center
  const double chordPoint = Point3::distance(point, center);
  const double chordCircle = Point3::distance(a, center);
  return chordPoint > (1.0 + 1e-6) * chordCircle;  // allow small sloppiness
}

bool checkFirstTriangulationIsDelaunay(const Point3 & p0, const Point3 & p1,
                                       const Point3 & p2, const Point3 & p3) {
  // First triangulation has triangles (p0,p1,p2) and (p2,p3,p0)
  // The Delaunay condition is that no other nodes lie in a triangle's circumcircle, so we check
  // that p3 is NOT in circumcircle of (p0,p1,p2) AND that p1 is NOT in circumcircle of (p2,p3,p0).
  return pointOutsideCircumcircle(p3, p0, p1, p2) && pointOutsideCircumcircle(p1, p2, p3, p0);
}
}  // namespace detail

namespace oops {

// -----------------------------------------------------------------------------

GeometryData::GeometryData(const atlas::FunctionSpace & fspace, const atlas::FieldSet & fset,
                           const bool topdown, const eckit::mpi::Comm & comm):
  fspace_(fspace), fset_(fset), comm_(comm), topdown_(topdown),
  firstTriangulationOfQuadsIsDelaunay_(),
  earth_(atlas::util::Earth::radius()), globalNodeTree_(earth_), localCellCenterTree_(earth_)
{
  // Exit constructor early if receiving an uninitialized or mesh-less FunctionSpace, because
  // the mesh and trees for the oops::UnstructuredInterpolator can't be built in this case.
  if (!fspace_) {
    Log::info() << "GeometryData::GeometryData received uninitialized FunctionSpace,"
      << " so skipping set up of interpolation data structures." << std::endl;
    return;
  } else if (fspace_.type() == "PointCloud" || fspace_.type() == "Spectral") {
    Log::info() << "GeometryData::GeometryData received FunctionSpace " << fspace_.type()
      << ", so skipping set up of interpolation data structures." << std::endl;
    return;
  }

  setMeshAndTriangulation();
  setLocalTree();
  setGlobalTree();
}

// -----------------------------------------------------------------------------

int GeometryData::closestTask(const double lat, const double lon) const {
  ASSERT(!globalNodeTree_.empty());
  atlas::PointLonLat target(lon, lat);
  target.normalise();
  const int itask = globalNodeTree_.closestPoint(target).payload();
  ASSERT(itask >= 0 && (size_t)itask < comm_.size());
  return itask;
}

// -----------------------------------------------------------------------------

bool GeometryData::containingTriangleAndBarycentricCoords(const double lat, const double lon,
    std::array<int, 3> & indices, std::array<double, 3> & baryCoords) const {
  ASSERT(mesh_);
  ASSERT(mesh_.nodes().has_field("xyz"));
  ASSERT(!localCellCenterTree_.empty());
  util::Timer timer("oops::GeometryData", "containingTriangleAndBarycentricCoords");

  const auto & connectivity = mesh_.cells().node_connectivity();
  const auto & xyz = atlas::array::make_view<double, 2>(mesh_.nodes().field("xyz"));

  const auto makePoint3 = [&](const int cell, const int node) -> atlas::Point3 {
    const int localindex = connectivity(cell, node);
    return atlas::Point3(xyz(localindex, 0), xyz(localindex, 1), xyz(localindex, 2));
  };

  const auto checkPointInSphericalTriangle = [&](const atlas::Point3 & point,
      const int cell, const int nodeA, const int nodeB, const int nodeC) -> bool {
    const atlas::Point3 a = makePoint3(cell, nodeA);
    const atlas::Point3 b = makePoint3(cell, nodeB);
    const atlas::Point3 c = makePoint3(cell, nodeC);
    const atlas::interpolation::element::Triag3D tri(a, b, c);
    const double sqrtArea = sqrt(tri.area());

    // Protect against the degenerate case where points coincide and the triangle has zero area.
    // This can occur for grids with singularities. If the triangle has zero area, we return early
    // and make sure it does not contain the target point.
    // Note: atlas computes XYZ mesh points using Earth geometry, so the 3D triangle lives on the
    // Earth's surface. We take the sqrt to go back to linear units, and compare to Earth's radius.
    if (sqrtArea < 1e-14 * earth_.radius()) {
      return false;
    }

    const atlas::interpolation::method::Ray ray(point);
    const double edgeEpsilon = 1e-15 * sqrtArea;
    const auto intersect = tri.intersects(ray, edgeEpsilon);
    if (intersect) {
      indices[0] = connectivity(cell, nodeA);
      indices[1] = connectivity(cell, nodeB);
      indices[2] = connectivity(cell, nodeC);
      baryCoords[0] = 1.0 - intersect.u - intersect.v;
      baryCoords[1] = intersect.u;
      baryCoords[2] = intersect.v;
      // The atlas coordinates u,v are in [0,1], but can still have roundoff-level negative 1-u-v
      if (baryCoords[0] < 0.0) {
        ASSERT(fabs(baryCoords[0]) < 1e-14);  // negative but larger than roundoff is a bug
        baryCoords[0] = 0.0;
      }
      return true;
    }
    return false;
  };

  // Find cell that contains target point
  atlas::PointLonLat pll(lon, lat);
  pll.normalise();
  atlas::Point3 p;
  earth_.lonlat2xyz(pll, p);

  // We check the 8 closest cells, following the example of the atlas unstructured interpolator
  const int nb_cells_to_check = std::min(8, mesh_.cells().size());
  const auto list = localCellCenterTree_.closestPoints(p, nb_cells_to_check);

  bool success = false;

  for (const auto & item : list) {
    const int cell = item.payload();
    const int nb_cols = connectivity.cols(cell);

    if (nb_cols == 3) {
      success = checkPointInSphericalTriangle(p, cell, 0, 1, 2);
      if (success) { break; }
    } else {
      ASSERT(!firstTriangulationOfQuadsIsDelaunay_.empty());
      if (firstTriangulationOfQuadsIsDelaunay_[cell]) {
        // triangle (p0,p1,p2)
        success = checkPointInSphericalTriangle(p, cell, 0, 1, 2);
        if (success) { break; }
        // triangle (p2,p3,p0)
        success = checkPointInSphericalTriangle(p, cell, 2, 3, 0);
        if (success) { break; }
      } else {
        // triangle (p3,p0,p1)
        success = checkPointInSphericalTriangle(p, cell, 3, 0, 1);
        if (success) { break; }
        // triangle (p1,p2,p3)
        success = checkPointInSphericalTriangle(p, cell, 1, 2, 3);
        if (success) { break; }
      }
    }
  }

  // The most likely explanation for failing to locate the target point is that it lies outside of
  // a regional grid. It is also possible (but unlikely) that the triangle containing the target
  // point is not in the group of triangles checked, i.e., is not one of the 8 closest triangles.
  // For now, we return false and let the client handle the failure to local the target point.
  return success;
}

// -----------------------------------------------------------------------------

void GeometryData::setGlobalTree() {
  ASSERT(globalNodeTree_.empty());
  util::Timer timer("oops::GeometryData", "setGlobalTree");

  // Count number of owned points
  const auto lonlat_view = atlas::array::make_view<double, 2>(fspace_.lonlat());
  const auto ghost_view = atlas::array::make_view<int, 1>(fspace_.ghost());
  const size_t nb_owned = [&ghost_view]() {
    int result = 0;
    for (atlas::idx_t jj = 0; jj < ghost_view.shape(0); ++jj) {
      if (ghost_view(jj) == 0) {
        ++result;
      }
    }
    return result;
  }();

  // Copy owned points into local buffer
  std::vector<double> lonlat(2 * nb_owned);
  size_t counter = 0;
  for (atlas::idx_t jj = 0; jj < ghost_view.shape(0); ++jj) {
    if (ghost_view(jj) == 0) {
      lonlat[2 * counter] = lonlat_view(jj, 0);
      lonlat[2 * counter + 1] = lonlat_view(jj, 1);
      ++counter;
    }
  }
  ASSERT(counter == nb_owned);

  // Collect global grid lats and lons
  const size_t nb_tasks = comm_.size();
  std::vector<size_t> sizes(nb_tasks);
  comm_.allGather(nb_owned, sizes.begin(), sizes.end());

  size_t nb_global = 0;
  for (size_t jtask = 0; jtask < nb_tasks; ++jtask) {
    nb_global += sizes[jtask];
  }

  std::vector<double> lonlat_global(2 * nb_global);
  mpi::allGatherv(comm_, lonlat, lonlat_global);

  // Arrange coordinates and task index for kd-tree
  std::vector<atlas::PointLonLat> nodes(nb_global);
  std::vector<int> tasks(nb_global);
  counter = 0;
  for (size_t jtask = 0; jtask < nb_tasks; ++jtask) {
    for (size_t jj = 0; jj < sizes[jtask]; ++jj) {
      nodes[counter] = atlas::PointLonLat(lonlat_global[2 * counter],
                                          lonlat_global[2 * counter + 1]);
      tasks[counter] = jtask;
      ++counter;
    }
  }
  ASSERT(counter == nb_global);

  // Create global kd-tree
  globalNodeTree_.build(nodes, tasks);
}

// -----------------------------------------------------------------------------

void GeometryData::setMeshAndTriangulation() {
  ASSERT(!mesh_);
  ASSERT(firstTriangulationOfQuadsIsDelaunay_.empty());
  util::Timer timer("oops::GeometryData", "setMeshAndTriangulation");

  // Setup mesh
  if (fspace_.type() == "NodeColumns") {
    const atlas::functionspace::NodeColumns nodecolumns(fspace_);
    mesh_ = nodecolumns.mesh();
  } else if (fspace_.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns structuredcolumns(fspace_);
    const atlas::StructuredGrid & grid = structuredcolumns.grid();
    if (fspace_.distribution() == "custom") {
      // Gather global partition field on root processor
      atlas::Field globalPartition = fspace_.createField<int>(
        atlas::option::name("partition") | atlas::option::global());
      fspace_.gather(fspace_.partition(), globalPartition);

      // Transform to a global partition vector
      std::vector<int> partition(grid.size());
      if (comm_.rank() == 0) {
        ASSERT(grid.size() == static_cast<int>(globalPartition.size()));
        const auto globalPartitionView = atlas::array::make_view<int, 1>(globalPartition);
        for (atlas::idx_t jj = 0; jj < grid.size(); ++jj) {
          partition[jj] = globalPartitionView(jj);
        }
      }

      // Broadcast global partition vector
      comm_.broadcast(partition, 0);

      // Create distribution and mesh
      atlas::grid::Distribution distribution;
      util::setupStructuredMeshWithCustomPartition(comm_, grid, partition, distribution, mesh_);
    } else {
      // Create mesh
      const atlas::StructuredMeshGenerator gen(grid.meshgenerator());
      mesh_ = gen(grid, atlas::grid::Partitioner(fspace_.distribution()));
    }

    // At this point, we have a mesh from the StructuredMeshGenerator that doesn't include a halo.
    // One can add a halo via action::build_mesh, but BEWARE several critical caveats:
    // 1. The halo added by build_halo is structured DIFFERENTLY than the halo in the
    //    StructuredColumns FunctionSpace. In other words: `fspace_.lonlat()` will be a different
    //    set of points from `mesh.nodes().lonlat()` -- the same owned points in the same order,
    //    but (in general) a different set of ghost points in a different order. Therefore, one
    //    CANNOT use a mesh-based computation to determine an index into the FunctionSpace/FieldSet,
    //    without first creating a mapping between the two halos.
    //    See https://github.com/JCSDA-internal/oops/issues/2621
    // 2. If the StructuredColumns FunctionSpace is global, then the build_halo action will NOT
    //    produce the expected halo surrounding every MPI partition -- in particular, it will not
    //    produce halos bridging the lon=0 meridian. It is not clear (as of atlas 0.37) if this is
    //    a bug or a feature, but the upshot is that for a global StructuredColumns the mesh halo
    //    is not complete and therefore likely unusable.
    //    See https://github.com/ecmwf/atlas/issues/200
    // To reduce the risk of incorrectly using the halos in the case of a structured mesh, we keep
    // the mesh halo-free for now.
  } else {
    ABORT(fspace_.type() + " function space not supported yet");
  }

  // Add XYZ field to mesh, because these 3d coordinates are used in interpolation setup
  atlas::mesh::actions::BuildXYZField()(mesh_);

  const size_t nb_cells = mesh_.cells().size();
  const auto & connectivity = mesh_.cells().node_connectivity();
  const auto & xyz = atlas::array::make_view<double, 2>(mesh_.nodes().field("xyz"));

  const auto makePoint3 = [&](const int cell, const int node) -> atlas::Point3 {
    const int localindex = connectivity(cell, node);
    return atlas::Point3(xyz(localindex, 0), xyz(localindex, 1), xyz(localindex, 2));
  };

  // If model has quad cells, triangulate them.
  //
  // We do this split to satisfy the Delaunay condition locally on the quad, giving the better of
  // the two possible splitting diagonals. Note this may not produce a globally-optimal Delaunay
  // triangulation, depending on how the model quads are set up.
  //
  // first triangulation  => split into triangles (p0,p1,p2) and (p2,p3,p0)
  // second triangulation => split into triangles (p3,p0,p1) and (p1,p2,p3)
  for (size_t i = 0; i < nb_cells; ++i) {
    const int nb_cols = connectivity.cols(i);
    ASSERT(nb_cols == 3 || nb_cols == 4);

    if (nb_cols == 4) {
      // Allocate triangulation data on first need
      if (firstTriangulationOfQuadsIsDelaunay_.empty()) {
        firstTriangulationOfQuadsIsDelaunay_.assign(nb_cells, false);
      }

      const atlas::Point3 p0 = makePoint3(i, 0);
      const atlas::Point3 p1 = makePoint3(i, 1);
      const atlas::Point3 p2 = makePoint3(i, 2);
      const atlas::Point3 p3 = makePoint3(i, 3);
      firstTriangulationOfQuadsIsDelaunay_[i] =
          detail::checkFirstTriangulationIsDelaunay(p0, p1, p2, p3);
    }
  }
}

// -----------------------------------------------------------------------------

void GeometryData::setLocalTree() {
  ASSERT(mesh_);
  ASSERT(mesh_.nodes().has_field("xyz"));
  ASSERT(localCellCenterTree_.empty());
  util::Timer timer("oops::GeometryData", "setLocalTree");

  const size_t nb_cells = mesh_.cells().size();

  const atlas::Field & centersField = atlas::mesh::actions::BuildCellCentres()(mesh_);
  const auto centersView = atlas::array::make_view<double, 2>(centersField);
  std::vector<atlas::Point3> centers(nb_cells);
  std::vector<int> indices(nb_cells);
  for (size_t i = 0; i < nb_cells; ++i) {
    centers[i] = atlas::Point3(centersView(i, 0), centersView(i, 1), centersView(i, 2));
    indices[i] = i;
  }
  localCellCenterTree_.build(centers, indices);
}

// -----------------------------------------------------------------------------

}  // namespace oops
