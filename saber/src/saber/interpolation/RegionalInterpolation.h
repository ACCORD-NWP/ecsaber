/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"
#include "atlas/redistribution/Redistribution.h"

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

#include "saber/interpolation/InterpElement.h"

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

class RegionalInterpolation {
 public:
  static const std::string classname()
    {return "quench::interpolation::RegionalInterpolation";}

  RegionalInterpolation(const atlas::FunctionSpace &,
                        const atlas::FunctionSpace &);
  ~RegionalInterpolation() {}

  void execute(const atlas::FieldSet &,
               atlas::FieldSet &) const;
  void executeAdjoint(atlas::FieldSet &,
                      const atlas::FieldSet &) const;

 private:
  void execute(const atlas::Field &,
               atlas::Field &) const;
  void executeAdjoint(atlas::Field &,
                      const atlas::Field &) const;

  const eckit::mpi::Comm & comm_;
  size_t srcSize_;
  std::vector<int> mpiTask_;
  size_t dstSize_;
  size_t srcSendSize_;
  size_t dstRecvSize_;
  std::vector<int> srcSendCounts_;
  std::vector<int> srcSendDispls_;
  std::vector<int> dstRecvCounts_;
  std::vector<int> dstRecvDispls_;
  std::vector<size_t> srcSendMapping_;
  std::vector<InterpElement> horInterp_;
};

}  // namespace interpolation
}  // namespace saber
