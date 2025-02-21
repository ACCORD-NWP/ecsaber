/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <Eigen/Dense>
#include <numeric>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

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

eckit::LocalConfiguration setObsValue(const std::string &);

// -----------------------------------------------------------------------------

eckit::LocalConfiguration setObsValue(const eckit::Configuration &,
                                      const std::string &);

// -----------------------------------------------------------------------------
// Timestamp
// -----------------------------------------------------------------------------

double timeStamp();

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
// Increment
// -----------------------------------------------------------------------------

template <typename MODEL>
void ones(oops::Increment<MODEL> & dx) {
  dx.increment().ones();
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
void maskAndSerialize(const oops::ObsVector<MODEL> & obsVec,
                      const oops::ObsVector<MODEL> & mask,
                      std::vector<double> & values) {
  obsVec.obsvector().maskAndSerialize(mask.obsvector(), values);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t serialSize(const oops::ObsVector<MODEL> & obsVec) {
  size_t len = obsVec.obsvector().serialSize();
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
void maskAndSerialize(const oops::ControlObsVector<MODEL> & ctlObsVec,
                      const oops::ControlObsVector<MODEL> & mask,
                      std::vector<double> & values) {
  util::maskAndSerialize(ctlObsVec.obsvector(), mask.obsvector(), values);
  if (ctlObsVec.hasBias()) {
    util::maskAndSerialize(ctlObsVec.biasvector(), mask.biasvector(), values);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t serialSize(const oops::ControlObsVector<MODEL> & ctlObsVec) {
  size_t len = util::serialSize(ctlObsVec.obsvector());
  if (ctlObsVec.hasBias()) {
    len += util::serialSize(ctlObsVec.biasvector());
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
  size_t len = 0;
  for (size_t idep = 0; idep < dep.size(); ++idep) {
    len += util::serialSize(dep[idep]);
  }

  std::vector<double> valid_values;
  valid_values.reserve(len);
  for (size_t idep = 0; idep < dep.size(); ++idep) {
    util::maskAndSerialize(dep[idep], mask[idep], valid_values);
  }
  const Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(valid_values.data(),
                                                          valid_values.size());
  return vec;
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
