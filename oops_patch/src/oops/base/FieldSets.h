/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/base/DataSetBase.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/Increment4D.h"
#include "oops/interface/Geometry.h"
#include "oops/util/ParallelFieldSetIO.h"
#include "oops/util/ECUtilities.h"

namespace oops {

// -----------------------------------------------------------------------------
class FieldSets : public DataSetBase<FieldSet3D, atlas::FunctionSpace> {
  typedef DataSetBase<FieldSet3D, atlas::FunctionSpace> Base_;

 public:
  /// @brief Creates a FieldSets from the IncrementSet. On creation fieldsets are
  ///        shared between Increment::fieldSet() and fieldsets in FieldSet4D.
  template<typename MODEL> FieldSets(const Geometry<MODEL> &,
                                     const JediVariables &,
                                     const std::vector<util::DateTime> &,
                                     const std::vector<int> & members,
                                     const std::vector<eckit::LocalConfiguration> &,
                                     const bool &);

  FieldSets(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
            const std::vector<int> &, const eckit::mpi::Comm &);

  FieldSets(const atlas::FunctionSpace &,
            const JediVariables &,
            const std::vector<util::DateTime> &,
            const eckit::Configuration &,
            const eckit::mpi::Comm &,
            const eckit::mpi::Comm & = oops::mpi::myself(),
            const eckit::mpi::Comm & = oops::mpi::myself());

  FieldSets(const atlas::FunctionSpace &,
            const JediVariables &,
            const util::ParallelFieldSetIO &,
            const std::vector<util::DateTime> &,
            const eckit::Configuration &,
            const eckit::mpi::Comm &,
            const eckit::mpi::Comm & = oops::mpi::myself(),
            const eckit::mpi::Comm & = oops::mpi::myself());

  /// @brief Emplace back FieldSet3D in empty FieldSets
  void emplace_back(const size_t &,
                    const size_t &,
                    const FieldSet3D &);

  /// @brief Multiplies each FieldSet3D in this FieldSets with the \p other.
  FieldSets & operator*=(const oops::FieldSet3D & other);
  FieldSets & operator*=(const double zz);

 private:
  std::string classname() const {return "FieldSets";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSets::FieldSets(const Geometry<MODEL> & geom,
                     const JediVariables & vars,
                     const std::vector<util::DateTime> & times,
                     const std::vector<int> & members,
                     const std::vector<eckit::LocalConfiguration> & memConfs,
                     const bool & removeMean)
  : Base_(times, eckit::mpi::self(), members, eckit::mpi::self()) {
  // Allocate ensemble
  for (size_t jm = 0; jm < members.size(); ++jm) {
    for (size_t jt = 0; jt < times.size(); ++jt) {
      this->dataset().emplace_back(std::make_unique<FieldSet3D>(times[jt], geom.geometry().getComm()));
    }
  }

  // Create mean pointer
  std::unique_ptr<oops::Increment4D<MODEL>> mean;

  for (size_t jm = 0; jm < members.size(); ++jm) {
    // Create 4D increment
    oops::Increment4D<MODEL> dx(geom, util::templatedVars<MODEL>(vars), times);

    // Read 4D increment
    if (memConfs[jm].has("states")) {
      eckit::LocalConfiguration memConf;
      memConf.set("increment", memConfs[jm].getSubConfiguration("states"));
      dx.read(memConf);
    } else {
      ASSERT(times.size() == 1);
      dx[0].read(memConfs[jm]);
    }

    if (removeMean) {
      if (!mean) {
        mean = std::make_unique<oops::Increment4D<MODEL>>(dx);
      } else {
        *mean += dx;
      }
    }

    // Add member
    for (size_t jt = 0; jt < times.size(); ++jt) {
      (*this)(jt, jm).shallowCopy(dx[jt].increment().fieldSet());
    }
  }

  if (removeMean) {
    // Normalize mean
    const double fact = 1.0 / static_cast<double>(this->ens_size());
    *mean *= fact;

    // Remove mean
    for (size_t jt = 0; jt < times.size(); ++jt) {
      FieldSet3D fsetMean(times[jt], geom.geometry().getComm());
      fsetMean.shallowCopy((*mean)[jt].increment().fieldSet());
      for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
        (*this)(jt, jm) -= fsetMean;
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
