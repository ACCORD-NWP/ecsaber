/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Interpolation.h"

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/util/FieldSetHelpers.h"

// -----------------------------------------------------------------------------

namespace quench {

// -----------------------------------------------------------------------------

Interpolation::Interpolation(const eckit::mpi::Comm & comm,
                             const atlas::grid::Partitioner & srcPartitioner,
                             const atlas::FunctionSpace & srcFspace,
                             const atlas::Grid & dstGrid,
                             const atlas::FunctionSpace & dstFspace)
  : comm_(comm), srcUid_(util::getGridUid(srcFspace)), dstUid_(dstGrid.uid()),
  dstFspace_(dstFspace), regionalInterpolation_(), globalInterpolation_() {
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;

  // Check domain type
  regionalSrcGrid_ = false;
  if (srcFspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns srcFs(srcFspace);
    regionalSrcGrid_ = !srcFs.grid().domain().global();
  }
  regionalDstGrid_ = false;
  if (dstFspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns dstFs(dstFspace);
    regionalDstGrid_ = !dstFs.grid().domain().global();
  }

  // Setup interpolation
  if (regionalSrcGrid_) {
    regionalInterpolation_ = std::make_shared<interp::RegionalInterpolation>(srcFspace, dstFspace);
  } else {
    // Configuration
    eckit::LocalConfiguration conf;
    if (regionalDstGrid_) {
      // Use OOPS unstructured interpolator
      conf.set("local interpolator type", "oops unstructured grid interpolator");
    } else {
      // Use ATLAS interpolator
      conf.set("local interpolator type", "atlas interpolator");
      eckit::LocalConfiguration interpMethodConf;
      interpMethodConf.set("type", "structured-linear2D");
      interpMethodConf.set("adjoint", true);
      conf.set("interpolation method", interpMethodConf);
    }

    // Empty fieldset
    atlas::FieldSet fset;

    // Geometry data
    oops::GeometryData geomData(srcFspace, fset, true, comm_);

    // Interpolator
    globalInterpolation_ = std::make_shared<oops::GlobalInterpolator>(conf, geomData, dstFspace,
      comm_);
  }

  // Test interpolation accuracy
  if (false) {
    // Create fields
    atlas::Field srcField = srcFspace.createField<double>(
      atlas::option::name("src") | atlas::option::levels(1));
    atlas::Field dstField = dstFspace.createField<double>(
      atlas::option::name("dst") | atlas::option::levels(1));

    // Define source field
    auto srcView = atlas::array::make_view<double, 2>(srcField);
    const auto srcLonLatView = atlas::array::make_view<double, 2>(srcFspace.lonlat());
    for (int jnode = 0; jnode < srcFspace.size(); ++jnode) {
      const double lon = srcLonLatView(jnode, 0);
      const double lat = srcLonLatView(jnode, 1);
      srcView(jnode, 0) = 0.5*(std::sin(2.0*M_PI*lon)*std::sin(2.0*M_PI*lat)+1.0);
    }

    // Create fieldsets
    atlas::FieldSet srcFieldSet;
    srcFieldSet.add(srcField);
    atlas::FieldSet dstFieldSet;
    dstFieldSet.add(dstField);

    // Interpolate
    execute(srcFieldSet, dstFieldSet);

    // Check accuracy
    double accuracy = 0.0;
    double maxVal = 0.0;
    double maxRefVal = 0.0;
    std::vector<double> locMax(2, 0.0);
    const auto dstView = atlas::array::make_view<double, 2>(dstField);
    const auto dstLonLatView = atlas::array::make_view<double, 2>(dstFspace.lonlat());
    for (int jnode = 0; jnode < dstFspace.size(); ++jnode) {
      const double lon = dstLonLatView(jnode, 0);
      const double lat = dstLonLatView(jnode, 1);
      const double refVal = 0.5*(std::sin(2.0*M_PI*lon)*std::sin(2.0*M_PI*lat)+1.0);
      const double diff = std::abs(dstView(jnode, 0)-refVal);
      if (diff > accuracy) {
        accuracy = diff;
        maxVal = dstView(jnode, 0);
        maxRefVal = refVal;
        locMax = {lon, lat};
      }
    }
    std::vector<double> accuracyPerTask(comm_.size());
    comm_.allGather(accuracy, accuracyPerTask.begin(), accuracyPerTask.end());
    size_t maxTask = std::distance(accuracyPerTask.begin(), std::max_element(
      accuracyPerTask.begin(), accuracyPerTask.end()));
    comm_.broadcast(maxVal, maxTask);
    comm_.broadcast(maxRefVal, maxTask);
    comm_.broadcast(locMax, maxTask);
    oops::Log::info() << std::setprecision(6) << "Info     :     Interpolation test accuracy: "
      << accuracy << " at " << locMax << " : " << std::setprecision(12) << maxVal << " != "
      << maxRefVal << std::endl;
  }

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::execute(const atlas::FieldSet & srcFieldSet,
                            atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  if (regionalSrcGrid_) {
    regionalInterpolation_->execute(srcFieldSet, targetFieldSet);
  } else {
    globalInterpolation_->apply(srcFieldSet, targetFieldSet);
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                   const atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  if (regionalSrcGrid_) {
    regionalInterpolation_->executeAdjoint(srcFieldSet, targetFieldSet);
  } else {
    globalInterpolation_->applyAD(srcFieldSet, targetFieldSet);
  }

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::insertVerticalInterpolation(const std::string & var,
                                                const std::vector<interp::InterpElement> & item) {
  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;

  if (verInterps_.find(var) != verInterps_.end()) {
    throw eckit::Exception("vertical interpolation already computed for this variables");
  }
  verInterps_.insert({var, item});

  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
