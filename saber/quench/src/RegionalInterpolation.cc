/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/RegionalInterpolation.h"

#include <iomanip>
#include <limits>

#include "atlas/array.h"
#include "atlas/grid/Distribution.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

namespace quench {

// -----------------------------------------------------------------------------

RegionalInterpolation::RegionalInterpolation(const atlas::FunctionSpace & srcFspace,
                                             const atlas::FunctionSpace & dstFspace)
  : comm_(eckit::mpi::comm()) {
  oops::Log::trace() << classname() << "::RegionalInterpolation starting" << std::endl;

  // Check source function space
  ASSERT(srcFspace.type() == "StructuredColumns");

  // Get grid parameters
  const atlas::functionspace::StructuredColumns srcFs(srcFspace);
  const atlas::RegularGrid srcGrid(srcFs.grid());
  const atlas::Projection & srcProj = srcGrid.projection();
  const size_t srcNx = srcGrid.nx();
  const size_t srcNy = srcGrid.ny();
  const double srcDx = srcGrid.dx();
  const double srcDy = std::abs(srcGrid.y(1)-srcGrid.y(0));
  const bool reversedY = srcGrid.y(1) < srcGrid.y(0);

  // Check grid regularity in y direction
  for (size_t srcJ = 0; srcJ < srcNy-1; ++srcJ) {
    if (reversedY) {
      ASSERT(std::abs(srcGrid.y(srcJ)-srcGrid.y(srcJ+1)-srcDy) < 1.0e-12*srcDy);
    } else {
      ASSERT(std::abs(srcGrid.y(srcJ+1)-srcGrid.y(srcJ)-srcDy) < 1.0e-12*srcDy);
    }
  }

  // Source grid indices
  const atlas::Field srcFieldIndexI = srcFs.index_i();
  const atlas::Field srcFieldIndexJ = srcFs.index_j();
  const auto srcIndexIView = atlas::array::make_view<int, 1>(srcFieldIndexI);
  const auto srcIndexJView = atlas::array::make_view<int, 1>(srcFieldIndexJ);
  srcSize_ = srcFs.size();

  // Destination grid size
  dstSize_ = dstFspace.size();

  // Ghost points
  const auto srcGhostView = atlas::array::make_view<int, 1>(srcFs.ghost());
  const auto dstGhostView = atlas::array::make_view<int, 1>(dstFspace.ghost());

  // Define reduced grid horizontal distribution
  std::vector<int> mpiTask(srcNx*srcNy, 0);
  for (size_t srcJnode = 0; srcJnode < srcSize_; ++srcJnode) {
    if (srcGhostView(srcJnode) == 0) {
      mpiTask[(srcIndexIView(srcJnode)-1)*srcNy+srcIndexJView(srcJnode)-1] = comm_.rank();
    }
  }
  comm_.allReduceInPlace(mpiTask.begin(), mpiTask.end(), eckit::mpi::sum());

  // Define local tree on destination grid
  std::vector<atlas::Point3> dstPoints;
  std::vector<size_t> dstIndices;
  const auto dstLonLatView = atlas::array::make_view<double, 2>(dstFspace.lonlat());
  for (size_t dstJnode = 0; dstJnode < dstSize_; ++dstJnode) {
    atlas::PointLonLat p({dstLonLatView(dstJnode, 0), dstLonLatView(dstJnode, 1)});
    srcProj.lonlat2xy(p);
    dstPoints.push_back(atlas::Point3(p[0], p[1], 0.0));
    dstIndices.push_back(dstJnode);
  }
  atlas::util::IndexKDTree dstTree;
  if (dstSize_ > 0) {
    dstTree.build(dstPoints, dstIndices);
  }
  const double radius = std::sqrt(srcDx*srcDx+srcDy*srcDy);

  // Delta for colocation
  const double eps = 1.0e-8;

  // RecvCounts and received points list
  dstRecvCounts_.resize(comm_.size());
  std::fill(dstRecvCounts_.begin(), dstRecvCounts_.end(), 0);
  std::vector<int> dstRecvPointsList;
  for (size_t srcJ = 0; srcJ < srcNy; ++srcJ) {
    double yMin, yMax;
    if (reversedY) {
      yMin = srcJ < srcNy-1 ? srcGrid.y(srcJ+1)-eps : -std::numeric_limits<double>::max();
      yMax = srcJ > 0 ? srcGrid.y(srcJ-1)+eps : std::numeric_limits<double>::max();
    } else {
      yMin = srcJ > 0 ? srcGrid.y(srcJ-1)-eps : -std::numeric_limits<double>::max();
      yMax = srcJ < srcNy-1 ? srcGrid.y(srcJ+1)+eps : std::numeric_limits<double>::max();
    }
    for (size_t srcI = 0; srcI < srcNx; ++srcI) {
      const double xMin = srcI > 0 ? srcGrid.x(srcI-1)-eps : -std::numeric_limits<double>::max();
      const double xMax = srcI < srcNx-1 ? srcGrid.x(srcI+1)+eps :
        std::numeric_limits<double>::max();

      bool pointsNeeded = false;
      if (dstSize_ > 0) {
        const atlas::Point3 p(srcGrid.x(srcI), srcGrid.y(srcJ), 0.0);
        const auto list = dstTree.closestPointsWithinRadius(p, radius);
        for (const auto & item : list) {
          const atlas::PointXYZ dstPoint = item.point();
          const size_t dstJnode = item.payload();
          if (dstGhostView(dstJnode) == 0) {
            const bool inX = (xMin <= dstPoint[0] && dstPoint[0] <= xMax);
            const bool inY = (yMin <= dstPoint[1] && dstPoint[1] <= yMax);
            if (inX && inY) {
              pointsNeeded = true;
              break;
            }
          }
        }
      }
      if (pointsNeeded) {
        ++dstRecvCounts_[mpiTask[srcI*srcNy+srcJ]];
        dstRecvPointsList.push_back(srcI*srcNy+srcJ);
      }
    }
  }

  // Buffer size
  dstRecvSize_ = dstRecvPointsList.size();

  if (dstRecvSize_ > 0) {
    // RecvDispls
    dstRecvDispls_.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      dstRecvDispls_.push_back(dstRecvDispls_[jt]+dstRecvCounts_[jt]);
    }

    // Allgather RecvCounts
    eckit::mpi::Buffer<int> dstRecvCountsBuffer(comm_.size());
    comm_.allGatherv(dstRecvCounts_.begin(), dstRecvCounts_.end(), dstRecvCountsBuffer);
    std::vector<int> dstRecvCountsGlb_ = std::move(dstRecvCountsBuffer.buffer);

    // SendCounts
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      srcSendCounts_.push_back(dstRecvCountsGlb_[jt*comm_.size()+comm_.rank()]);
    }

    // Buffer size
    srcSendSize_ = 0;
    for (const auto & n : srcSendCounts_) srcSendSize_ += n;

    // SendDispls
    srcSendDispls_.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      srcSendDispls_.push_back(srcSendDispls_[jt]+srcSendCounts_[jt]);
    }

    // Ordered received points list
    std::vector<size_t> dstRecvOffset(comm_.size(), 0);
    std::vector<int> dstRecvPointsListOrdered(dstRecvSize_);
    for (size_t jr = 0; jr < dstRecvSize_; ++jr) {
      const size_t srcI = dstRecvPointsList[jr]/srcNy;
      const size_t srcJ = dstRecvPointsList[jr]-srcI*srcNy;
      size_t jt = mpiTask[srcI*srcNy+srcJ];
      size_t jro = dstRecvDispls_[jt]+dstRecvOffset[jt];
      dstRecvPointsListOrdered[jro] = dstRecvPointsList[jr];
      ++dstRecvOffset[jt];
    }
    std::vector<int> srcSentPointsList(srcSendSize_);
    comm_.allToAllv(dstRecvPointsListOrdered.data(), dstRecvCounts_.data(), dstRecvDispls_.data(),
                    srcSentPointsList.data(), srcSendCounts_.data(), srcSendDispls_.data());

    // Sort indices
    std::vector<int> gij;
    for (size_t srcJnode = 0; srcJnode < srcSize_; ++srcJnode) {
      if (srcGhostView(srcJnode) == 0) {
        gij.push_back((srcIndexIView(srcJnode)-1)*srcNy+srcIndexJView(srcJnode)-1);
      } else {
        gij.push_back(-1);
      }
    }
    std::vector<size_t> gidx(srcSize_);
    std::iota(gidx.begin(), gidx.end(), 0);
    std::stable_sort(gidx.begin(), gidx.end(), [&gij](size_t i1, size_t i2)
      {return gij[i1] < gij[i2];});
    std::vector<size_t> ridx(srcSendSize_);
    std::iota(ridx.begin(), ridx.end(), 0);
    std::stable_sort(ridx.begin(), ridx.end(), [&srcSentPointsList](size_t i1, size_t i2)
      {return srcSentPointsList[i1] < srcSentPointsList[i2];});

    // Mapping for sent points
    srcSendMapping_.resize(srcSendSize_);
    size_t srcJnode = 0;
    for (size_t js = 0; js < srcSendSize_; ++js) {
      while (gij[gidx[srcJnode]] < srcSentPointsList[ridx[js]]) {
        ++srcJnode;
        ASSERT(srcJnode < srcSize_);
      }
      srcSendMapping_[ridx[js]] = gidx[srcJnode];
    }

    // Sort indices
    std::vector<size_t> idx(dstRecvPointsListOrdered.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&dstRecvPointsListOrdered](size_t i1, size_t i2)
      {return dstRecvPointsListOrdered[i1] < dstRecvPointsListOrdered[i2];});

    // Compute horizontal interpolation
    for (size_t dstJnode = 0; dstJnode < dstSize_; ++dstJnode) {
      // Interpolation element default values
      std::vector<std::pair<size_t, double>> operations;

      if (dstGhostView(dstJnode) == 0) {
        // Destination grid indices
        const double dstX = dstPoints[dstJnode][0];
        bool colocatedX = false;
        int indexI = -1;
        for (size_t srcI = 0; srcI < srcNx-1; ++srcI) {
          if (std::abs(dstX-srcGrid.x(srcI)) < eps) {
            indexI = srcI;
            colocatedX = true;
          }
          if (srcGrid.x(srcI)+eps < dstX && dstX < srcGrid.x(srcI+1)-eps) {
            indexI = srcI;
            colocatedX = false;
          }
        }
        if (std::abs(dstX-srcGrid.x(srcNx-1)) < eps) {
          indexI = srcNx-1;
          colocatedX = true;
        }
        const double dstY = dstPoints[dstJnode][1];
        bool colocatedY = false;
        int indexJ = -1;
        for (size_t srcJ = 0; srcJ < srcNy-1; ++srcJ) {
          if (std::abs(dstY-srcGrid.y(srcJ)) < eps) {
            indexJ = srcJ;
            colocatedY = true;
          }
          if (reversedY) {
            if (srcGrid.y(srcJ+1)+eps < dstY && dstY < srcGrid.y(srcJ)-eps) {
              indexJ = srcJ;
              colocatedY = false;
            }
          } else {
            if (srcGrid.y(srcJ)+eps < dstY && dstY < srcGrid.y(srcJ+1)-eps) {
              indexJ = srcJ;
              colocatedY = false;
            }
          }
        }
        if (std::abs(dstY-srcGrid.y(srcNy-1)) < eps) {
          indexJ = srcNy-1;
          colocatedY = true;
        }

        if (indexI == -1 || indexJ == -1) {
          // Point outside of the domain, using nearest neighbor
          if (indexI > -1) {
            if (!colocatedX &&
              (std::abs(dstX-srcGrid.x(indexI+1)) < std::abs(dstX-srcGrid.x(indexI)))) {
              indexI += 1;
            }
          } else {
            if (std::abs(dstX-srcGrid.x(0)) < std::abs(dstX-srcGrid.x(srcNx-1))) {
              indexI = 0;
            } else {
              indexI = srcNx-1;
            }
          }
          if (indexJ > -1) {
            if (!colocatedY &&
              (std::abs(dstY-srcGrid.y(indexJ+1)) < std::abs(dstY-srcGrid.y(indexJ)))) {
              indexJ += 1;
            }
          } else {
            if (std::abs(dstY-srcGrid.y(0)) < std::abs(dstY-srcGrid.y(srcNy-1))) {
              indexJ = 0;
            } else {
              indexJ = srcNy-1;
            }
            std::cout << "WARNING: point outside of the domain" << std::endl;
          }

          // Colocated point (actually nearest neighbor)
          colocatedX = true;
          colocatedY = true;
        }

        // Bilinear interpolation factor
        const double alphaX = 1.0-(srcGrid.x(indexI)+srcDx-dstX)/srcDx;
        const double alphaY = reversedY ? (srcGrid.y(indexJ)-dstY)/srcDy
          : 1.0-(srcGrid.y(indexJ)+srcDy-dstY)/srcDy;

        // Points to find
        std::vector<bool> toFind = {true, !colocatedX, !colocatedY, !colocatedX && !colocatedY};
        std::vector<size_t> valueToFind = {indexI*srcNy+indexJ, (indexI+1)*srcNy+indexJ,
          indexI*srcNy+(indexJ+1), (indexI+1)*srcNy+(indexJ+1)};
        std::vector<int> foundIndex(4, -1);

        // Binary search for each point
        for (size_t jj = 0; jj < 4; ++jj) {
          if (toFind[jj]) {
            size_t low = 0;
            size_t high = dstRecvPointsListOrdered.size()-1;
            while (low <= high) {
              size_t mid = low+(high-low)/2;
              if (valueToFind[jj] == static_cast<size_t>(dstRecvPointsListOrdered[idx[mid]])) {
                foundIndex[jj] = idx[mid];
                break;
              }
              if (valueToFind[jj] > static_cast<size_t>(dstRecvPointsListOrdered[idx[mid]])) {
                low = mid+1;
              }
              if (valueToFind[jj] < static_cast<size_t>(dstRecvPointsListOrdered[idx[mid]])) {
                high = mid-1;
              }
            }
            ASSERT(foundIndex[jj] > -1);
            ASSERT(static_cast<size_t>(dstRecvPointsListOrdered[foundIndex[jj]]) ==
              valueToFind[jj]);
          }
        }

        // Create interpolation operations
        if (colocatedX && colocatedY) {
          // Colocated point
          operations.push_back(std::make_pair(foundIndex[0], 1.0));
        } else if (colocatedY) {
          // Linear interpolation along x
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaX));
          operations.push_back(std::make_pair(foundIndex[1], alphaX));
        } else if (colocatedX) {
          // Linear interpolation along y
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaY));
          operations.push_back(std::make_pair(foundIndex[2], alphaY));
        } else {
          // Bilinear interpolation
          operations.push_back(std::make_pair(foundIndex[0], (1.0-alphaX)*(1.0-alphaY)));
          operations.push_back(std::make_pair(foundIndex[1], alphaX*(1.0-alphaY)));
          operations.push_back(std::make_pair(foundIndex[2], (1.0-alphaX)*alphaY));
          operations.push_back(std::make_pair(foundIndex[3], alphaX*alphaY));
        }
      }
      horInterp_.push_back(InterpElement(operations));
    }
  }

  oops::Log::trace() << classname() << "::RegionalInterpolation done" << std::endl;
}

//-------------------------------------------------------------------------------------------------

void RegionalInterpolation::execute(const atlas::Field & srcField,
                                    atlas::Field & dstField) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  // Check number of levels
  ASSERT(srcField.levels() == dstField.levels());
  const size_t nz = srcField.levels();

  // Scale counts and displs for all levels
  std::vector<int> srcSendCounts3D(comm_.size());
  std::vector<int> srcSendDispls3D(comm_.size());
  std::vector<int> dstRecvCounts3D(comm_.size());
  std::vector<int> dstRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    srcSendCounts3D[jt] = srcSendCounts_[jt]*nz;
    srcSendDispls3D[jt] = srcSendDispls_[jt]*nz;
    dstRecvCounts3D[jt] = dstRecvCounts_[jt]*nz;
    dstRecvDispls3D[jt] = dstRecvDispls_[jt]*nz;
  }

  // Serialize
  const auto srcView = atlas::array::make_view<double, 2>(srcField);
  std::vector<double> srcSendVec(srcSendSize_*nz);
  for (size_t js = 0; js < srcSendSize_; ++js) {
    size_t srcJnode = srcSendMapping_[js];
    for (size_t k = 0; k < nz; ++k) {
      ASSERT(js*nz+k < srcSendSize_*nz);
      srcSendVec[js*nz+k] = srcView(srcJnode, k);
    }
  }

  // Communication
  std::vector<double> dstRecvVec(dstRecvSize_*nz);

  comm_.allToAllv(srcSendVec.data(), srcSendCounts3D.data(), srcSendDispls3D.data(),
                  dstRecvVec.data(), dstRecvCounts3D.data(), dstRecvDispls3D.data());

  // Interpolation
  const auto dstGhostView = atlas::array::make_view<int, 1>(dstField.functionspace().ghost());
  auto dstView = atlas::array::make_view<double, 2>(dstField);
  dstView.assign(0.0);
  for (size_t dstJnode = 0; dstJnode < dstSize_; ++dstJnode) {
    if (dstGhostView(dstJnode) == 0) {
      for (const auto & horOperation : horInterp_[dstJnode].operations()) {
        for (size_t k0 = 0; k0 < nz; ++k0) {
          size_t mIndex = horOperation.first*nz+k0;
          dstView(dstJnode, k0) += horOperation.second*dstRecvVec[mIndex];
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void RegionalInterpolation::executeAdjoint(atlas::Field & srcField,
                                           const atlas::Field & dstField) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  // Check number of levels
  ASSERT(srcField.levels() == dstField.levels());
  const size_t nz = srcField.levels();

  // Scale counts and displs for all levels
  std::vector<int> srcSendCounts3D(comm_.size());
  std::vector<int> srcSendDispls3D(comm_.size());
  std::vector<int> dstRecvCounts3D(comm_.size());
  std::vector<int> dstRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    srcSendCounts3D[jt] = srcSendCounts_[jt]*nz;
    srcSendDispls3D[jt] = srcSendDispls_[jt]*nz;
    dstRecvCounts3D[jt] = dstRecvCounts_[jt]*nz;
    dstRecvDispls3D[jt] = dstRecvDispls_[jt]*nz;
  }

  // Copy destination field
  atlas::Field dstTmpField = dstField.clone();

  // Halo exchange
  dstTmpField.adjointHaloExchange();

  // Interpolation adjoint
  const auto dstGhostView = atlas::array::make_view<int, 1>(dstField.functionspace().ghost());
  const auto dstView = atlas::array::make_view<double, 2>(dstTmpField);
  std::vector<double> dstRecvVec(dstRecvSize_*nz, 0.0);
  for (size_t dstJnode = 0; dstJnode < dstSize_; ++dstJnode) {
    if (dstGhostView(dstJnode) == 0) {
      for (const auto & horOperation : horInterp_[dstJnode].operations()) {
        for (size_t k0 = 0; k0 < nz; ++k0) {
          size_t mIndex = horOperation.first*nz+k0;
          dstRecvVec[mIndex] += horOperation.second*dstView(dstJnode, k0);
        }
      }
    }
  }

  // Communication
  std::vector<double> srcSendVec(srcSendSize_*nz);
  comm_.allToAllv(dstRecvVec.data(), dstRecvCounts3D.data(), dstRecvDispls3D.data(),
                  srcSendVec.data(), srcSendCounts3D.data(), srcSendDispls3D.data());

  // Deserialize
  auto srcView = atlas::array::make_view<double, 2>(srcField);
  srcView.assign(0.0);
  for (size_t js = 0; js < srcSendSize_; ++js) {
    size_t srcJnode = srcSendMapping_[js];
    for (size_t k = 0; k < nz; ++k) {
      srcView(srcJnode, k) += srcSendVec[js*nz+k];
    }
  }

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

void RegionalInterpolation::execute(const atlas::FieldSet & srcFieldSet,
                                    atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  srcFieldSet.haloExchange();
  for (auto & srcField : srcFieldSet) {
    execute(srcField, targetFieldSet[srcField.name()]);
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void RegionalInterpolation::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                           const atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  for (auto & srcField : srcFieldSet) {
    executeAdjoint(srcField, targetFieldSet[srcField.name()]);
  }
  srcFieldSet.adjointHaloExchange();
  srcFieldSet.set_dirty();

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
