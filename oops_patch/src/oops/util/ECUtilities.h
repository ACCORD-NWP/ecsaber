/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Dense>
#include <numeric>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlObsVector.h"
#include "oops/base/Departures.h"
#include "oops/base/Increment4D.h"
#include "oops/base/Observations.h"
#include "oops/base/Variables.h"
#include "oops/generic/ObsErrorDiag.h"
#include "oops/interface/ObsVector.h"

namespace util {

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration &,
               const size_t &);

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

// -----------------------------------------------------------------------------

void expandEnsembleTemplate(eckit::LocalConfiguration &,
                            const size_t &);

// -----------------------------------------------------------------------------
// Variables
// -----------------------------------------------------------------------------

eckit::LocalConfiguration templatedVarsConf(const oops::JediVariables &);

// -----------------------------------------------------------------------------
// Increment4D
// -----------------------------------------------------------------------------

template<typename MODEL>
void dirac4D(const eckit::Configuration & conf,
             oops::Increment4D<MODEL> & incr4D) {
  if (incr4D.first() == incr4D.last()) {
    incr4D[0].increment().dirac(conf);
  } else {
    const std::vector<eckit::LocalConfiguration> confs = conf.getSubConfigurations();
    ASSERT(incr4D.last()-incr4D.first()+1 == confs.size());
    for (int jt = incr4D.first(); jt <= incr4D.last(); ++jt) {
      if (!confs[jt].empty()) {
        incr4D[jt].increment().dirac(confs[jt]);
      } else {
        incr4D[jt].increment().zero();
      }
    }
  }
}

// -----------------------------------------------------------------------------
// ObsVector
// -----------------------------------------------------------------------------

template <typename MODEL>
void ones(oops::ObsVector<MODEL> & obsVec) {
  obsVec.obsvector().ones();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
unsigned int nobs(const oops::ObsVector<MODEL> & obsVec) {
  int nobs = obsVec.obsvector().nobs();
  return nobs;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void mask(oops::ObsVector<MODEL> & obsVec,
          const oops::ObsVector<MODEL> & mask) {
  obsVec.obsvector().mask(mask.obsvector());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Eigen::VectorXd packEigen(const oops::ObsVector<MODEL> & obsVec,
                          const oops::ObsVector<MODEL> & mask) {
  Eigen::VectorXd vec = obsVec.obsvector().packEigen(mask.obsvector());
  return vec;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t packEigenSize(const oops::ObsVector<MODEL> & obsVec,
                     const oops::ObsVector<MODEL> & mask) {
  size_t len = obsVec.obsvector().packEigenSize(mask.obsvector());
  return len;
}

// -----------------------------------------------------------------------------
// ControlObsVector
// -----------------------------------------------------------------------------

template <typename MODEL>
void ones(oops::ControlObsVector<MODEL> & ctlObsVec) {
  util::ones(ctlObsVec.obsvector());
  if (ctlObsVec.hasBias()) util::ones(ctlObsVec.biasvector());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void mask(oops::ControlObsVector<MODEL> & ctlObsVec,
          const oops::ControlObsVector<MODEL> & mask) {
  util::mask(ctlObsVec.obsvector(), mask.obsvector());
  if (ctlObsVec.hasBias()) {
    util::mask(ctlObsVec.biasvector(), mask.biasvector());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Eigen::VectorXd packEigen(const oops::ControlObsVector<MODEL> & ctlObsVec,
                          const oops::ControlObsVector<MODEL> & mask) {
  Eigen::VectorXd vec = util::packEigen(ctlObsVec.obsvector(), mask.obsvector());
  if (ctlObsVec.hasBias()) {
    Eigen::VectorXd vecBias = util::packEigen(ctlObsVec.biasvector(), mask.biasvector());
    Eigen::VectorXd vecJoined(vec.size() + vecBias.size());
    vecJoined << vec, vecBias;
    return vecJoined;
  } else {
    return vec;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t packEigenSize(const oops::ControlObsVector<MODEL> & ctlObsVec,
                     const oops::ControlObsVector<MODEL> & mask) {
  size_t len = util::packEigenSize(ctlObsVec.obsvector(), mask.obsvector());
  if (ctlObsVec.hasBias()) {
    len += util::packEigenSize(ctlObsVec.biasvector(), mask.biasvector());
  }
  return len;
}

// -----------------------------------------------------------------------------
// Departures
// -----------------------------------------------------------------------------

template <typename MODEL>
void ones(oops::Departures<MODEL> & dep) {
  for (std::size_t jj = 0; jj < dep.size(); ++jj) {
    util::ones(dep[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void mask(oops::Departures<MODEL> & dep,
          const oops::Departures<MODEL> & mask) {
  for (size_t ii = 0; ii < dep.size(); ++ii) {
    util::mask(dep[ii], mask[ii]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Eigen::VectorXd packEigen(const oops::Departures<MODEL> & dep,
                          const oops::Departures<MODEL> & mask) {
  std::vector<size_t> len(dep.size());
  for (size_t idep = 0; idep < dep.size(); ++idep) {
    len[idep] = util::packEigenSize(dep[idep], mask[idep]);
  }
  size_t all_len = std::accumulate(len.begin(), len.end(), 0);

  Eigen::VectorXd vec(all_len);
  size_t ii = 0;
  for (size_t idep = 0; idep < dep.size(); ++idep) {
    vec.segment(ii, len[idep]) = util::packEigen(dep[idep], mask[idep]);
    ii += len[idep];
  }
  return vec;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t packEigenSize(const oops::Departures<MODEL> & dep,
                     const oops::Departures<MODEL> & mask) {
  size_t len = 0;
  for (size_t idep = 0; idep < dep.size(); ++idep) {
    len += util::packEigenSize(dep[idep], mask[idep]);
  }
  return len;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double rms(oops::Departures<MODEL> & dep) {
  double zz = dep.dot_product_with(dep);
  unsigned int n = 0;
  for (std::size_t jj = 0; jj < dep.size(); ++jj) {
    n += dep[jj].size();
  }
  if (n > 0) zz = zz/static_cast<double>(n);
  zz = std::sqrt(zz);
  return zz;
}


// -----------------------------------------------------------------------------
// Observations
// -----------------------------------------------------------------------------

template <typename MODEL>
void ones(oops::Observations<MODEL> & obs) {
  for (std::size_t jj = 0; jj < obs.size(); ++jj) {
    obs[jj]->ones();
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
