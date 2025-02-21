/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENTENSEMBLE4D_H_
#define OOPS_BASE_INCREMENTENSEMBLE4D_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "oops/interface/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ECUtilities.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of 4D increments
template<typename MODEL> class IncrementEnsemble4D {
  typedef Geometry<MODEL>            Geometry_;
  typedef typename MODEL::GeometryIterator  GeometryIterator__;
  typedef State4D<MODEL>             State4D_;
  typedef StateEnsemble4D<MODEL>     StateEnsemble4D_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef Variables<MODEL>           Variables_;

 public:
  /// Constructor
  IncrementEnsemble4D(const Geometry_ & resol,
                      const JediVariables & vars,
                      const std::vector<util::DateTime> &,
                      const int rank);
  /// \brief construct ensemble of perturbations as \p ens - \p mean; holding
  //         \p vars variables
  IncrementEnsemble4D(const StateEnsemble4D_ & ens, const State4D_ & mean,
                      const JediVariables & vars);

  /// Accessors
  size_t size() const {return ensemblePerturbs_.size();}
  Increment4D_ & operator[](const size_t ii) {return ensemblePerturbs_[ii];}
  const Increment4D_ & operator[](const size_t ii) const {return ensemblePerturbs_[ii];}

  /// Eigen interface
  void packEigen(Eigen::MatrixXd &, const GeometryIterator__ &, const size_t &) const;
  void setEigen(const Eigen::MatrixXd &, const GeometryIterator__ &, const size_t &);

 private:
  std::vector<Increment4D_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsemble4D<MODEL>::IncrementEnsemble4D(const Geometry_ & resol, const JediVariables & vars,
                                                const std::vector<util::DateTime> & timeslots,
                                                const int rank)
  : ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(rank);
  const Variables_ varsT(templatedVarsConf(vars));
  for (int m = 0; m < rank; ++m) {
    ensemblePerturbs_.emplace_back(resol, varsT, timeslots);
  }
  Log::trace() << "IncrementEnsemble4D:contructor done" << std::endl;
}

// ====================================================================================

template<typename MODEL>
IncrementEnsemble4D<MODEL>::IncrementEnsemble4D(const StateEnsemble4D_ & ensemble,
                                                const State4D_ & mean, const JediVariables & vars)
  : ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(ensemble.size());
  const Variables_ varsT(templatedVarsConf(vars));
  for (size_t ii = 0; ii < ensemble.size(); ++ii) {
    ensemblePerturbs_.emplace_back(ensemble[ii][0].geometry(), varsT,
                                   ensemble[ii].times());
    ensemblePerturbs_[ii].diff(ensemble[ii], mean);
  }
  Log::trace() << "IncrementEnsemble4D:contructor(StateEnsemble4D) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble4D<MODEL>::packEigen(Eigen::MatrixXd & X,
                                         const GeometryIterator__ & gi,
                                         const size_t & itime) const
{
  size_t ngp = ensemblePerturbs_[0][itime].increment().getLocal(gi).getVals().size();
  size_t nens = ensemblePerturbs_.size();
  X.resize(ngp, nens);
  for (size_t iens=0; iens < nens; ++iens) {
    LocalIncrement gp = ensemblePerturbs_[iens][itime].increment().getLocal(gi);
    std::vector<double> tmp1 = gp.getVals();
    for (size_t iv=0; iv < ngp; ++iv) {
      X(iv, iens) = tmp1[iv];
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble4D<MODEL>::setEigen(const Eigen::MatrixXd & X,
                                        const GeometryIterator__ & gi,
                                        const size_t & itime)
{
  size_t ngp = ensemblePerturbs_[0][itime].increment().getLocal(gi).getVals().size();
  size_t nens = ensemblePerturbs_.size();

  LocalIncrement gptmp = ensemblePerturbs_[0][itime].increment().getLocal(gi);
  std::vector<double> tmp = gptmp.getVals();

  for (size_t iens=0; iens < nens; ++iens) {
    for (size_t iv=0; iv < ngp; ++iv) {
      tmp[iv] = X(iv, iens);
    }
    gptmp.setVals(tmp);
    ensemblePerturbs_[iens][itime].increment().setLocal(gptmp, gi);
  }
}

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE4D_H_
