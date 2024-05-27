/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iomanip>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/functionspace.h"

#include "oops/generic/GlobalInterpolator.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

#include "saber/interpolation/InterpElement.h"
#include "saber/interpolation/RegionalInterpolation.h"

namespace atlas {
  class Field;
  class Grid;
  namespace grid {
    class Partitioner;
  }
}

namespace interp = saber::interpolation;

namespace quench {

// -----------------------------------------------------------------------------

class Interpolation {
 public:
  static const std::string classname()
    {return "quench::interpolation::Interpolation";}

  Interpolation(const eckit::mpi::Comm &,
                const atlas::grid::Partitioner &,
                const atlas::FunctionSpace &,
                const atlas::Grid &,
                const atlas::FunctionSpace &);
  ~Interpolation() {}

  void execute(const atlas::FieldSet &,
               atlas::FieldSet &) const;
  void executeAdjoint(atlas::FieldSet &,
                      const atlas::FieldSet &) const;

  void insertVerticalInterpolation(const std::string &,
                                   const std::vector<interp::InterpElement> &);
  std::vector<interp::InterpElement> & verticalInterpolation(const std::string & var)
    {return verInterps_[var];}

  const std::string & srcUid() const
    {return srcUid_;}
  const std::string & dstUid() const
    {return dstUid_;}
  const atlas::FunctionSpace & dstFspace() const
    {return dstFspace_;}

 private:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Grids UID
  std::string srcUid_;
  std::string dstUid_;

  // Destination function space
  atlas::FunctionSpace dstFspace_;

  // Grid type
  bool regionalSrcGrid_;
  bool regionalDstGrid_;

  // Regional interpolation
  std::shared_ptr<interp::RegionalInterpolation> regionalInterpolation_;

  // Global interpolation
  std::shared_ptr<oops::GlobalInterpolator> globalInterpolation_;

  // Vertical interpolations
  std::unordered_map<std::string, std::vector<interp::InterpElement>> verInterps_;
};

}  // namespace quench
