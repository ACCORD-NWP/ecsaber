/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/assimilation/ControlVector.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/JbMatrix.h"
#include "oops/assimilation/SQRTMinimizer.h"
#include "oops/assimilation/SpectralSqrtLMP.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/assimilation/UtHtRinvHUMatrix.h"
#include "oops/util/Logger.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/formats.h"

#include "saber/oops/RitzPairs.h"

#include "util/dot_product.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SQRTBPLanczosEVILMinimizer : public oops::SQRTMinimizer<MODEL> {
  using CostFct_ = oops::CostFunction<MODEL>;
  using CtrlVec_ = oops::ControlVector<MODEL>;
  using Jbmat_ = oops::JbMatrix<MODEL>;
  using UtHtRinvHU_ = oops::UtHtRinvHUMatrix<MODEL>;

 public:
  const std::string classname() const final { return "SQRTBPLanczosEVILMinimizer"; }

  /// Constructor, destructor
  SQRTBPLanczosEVILMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~SQRTBPLanczosEVILMinimizer() {}

 private:
  double solve(CtrlVec_ &, CtrlVec_ &, const Jbmat_ &, const UtHtRinvHU_ &,
               const bool &, const size_t &, const size_t &, const double &,
               const double &, const double &) final;

  void releaseResources();

  /// Checkpointing, restart of the limited memory preconditioner
  void checkpointLMP(eckit::LocalConfiguration &) const final;
  void restartLMP(const eckit::Configuration &) final;

  /// Quadratic cost function calculations
  void setupQuadCost(const double &costJ0Jb, const double &costJ0JoJc);
  void calcQuadCost(const CtrlVec_ &, const std::vector<double> &);
  void printQuadCost(const size_t &) const;

  /// Other diagnostics
  bool calcRitzInformation();
  void printRitzInformation(const size_t &) const;

  /// Data members
  const eckit::LocalConfiguration conf_;
  const CostFct_ &J_;

  /// LMP
  oops::SpectralSqrtLMP<MODEL> lmp_;

  /// Local
  RitzPairs<CtrlVec_> ritzPairs_;
  std::vector<std::shared_ptr<CtrlVec_> > vvecs_;
  std::vector<std::shared_ptr<CtrlVec_> > zvecs_;
  std::vector<double> alphas_;
  std::vector<double> betas_;

  std::vector<double> evals_;
  std::vector<double> eigvals_;
  std::vector<double> erritz_;
  std::vector<double> erritzlm_;
  std::vector<std::vector<double> > evecs_;

  std::vector<std::string> error_messages_;

  /// Quadratic cost function calculations data members
  double costJ0_;
  double costJ0Jb_;
  double costJ0JoJc_;

  double costJ_;
  double costJb_;
  double costJbLocal_;  // component of costJb from current minimization
  double costJoJc_;
};

// =============================================================================

template <typename MODEL>
SQRTBPLanczosEVILMinimizer<MODEL>::SQRTBPLanczosEVILMinimizer(
    const eckit::Configuration &conf, const CostFct_ &J)
    : oops::SQRTMinimizer<MODEL>(conf, J),
      conf_(conf),
      J_(J),
      lmp_(conf, J),
      vvecs_(),
      alphas_(),
      betas_(),
      costJ0_(0),
      costJ0Jb_(0),
      costJ0JoJc_(0),
      costJ_(0),
      costJb_(0),
      costJbLocal_(0),
      costJoJc_(0) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
double SQRTBPLanczosEVILMinimizer<MODEL>::solve(
    CtrlVec_ &dv, CtrlVec_ &rr, const Jbmat_ &Jb, const UtHtRinvHU_ &UtHtRinvHU,
    const bool &linfoconv, const size_t &maxiter, const size_t &miniter,
    const double &tolerance, const double &costJ0Jb, const double &costJ0JoJc) {
  // dv   control vector

  CtrlVec_ dv0(dv);
  CtrlVec_ jbzz(dv, false);
  CtrlVec_ vv(rr, false);
  CtrlVec_ zz(rr, false);

  std::vector<double> ss;
  std::vector<double> dd;

  // J0
  this->setupQuadCost(costJ0Jb, costJ0JoJc);

  // Change resolution of LMP vectors
  lmp_.changeResolution();

  // Postprocess preconditioning vectors
  lmp_.postprocess();

  // P rr_{0}
  lmp_.inverseMultiply(rr, zz);

  // beta_{0} = sqrt( rr_{0}^T zz_{0} )
  double beta = sqrt(dot_product(rr, zz));
  const double beta0 = beta;

  // v_{1} = r_{0} / beta_{0}
  vv = rr;
  vv *= 1.0 / beta;

  // z_{1} = z_{0} / beta_{0}
  zz *= 1.0 / beta;

  // vvecs[0] = v_{1} ---> for re-orthogonalization
  vvecs_.push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(vv)));
  zvecs_.push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(zz)));

  ritzPairs_.vVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(vv)));
  ritzPairs_.zVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(zz)));

  double normReduction = 1.0;
  double zdfi = 0.0;
  bool convergenceReached = false;
  size_t jiter = 0;

  while (jiter < maxiter) {
    oops::Log::info() << "SQRTBPLanczos Starting Iteration " << jiter + 1
                      << std::endl;

    // vv_{i+1} = ( I + U^T H^T R^{-1} H U ) zz_{i}
    UtHtRinvHU.multiply(zz, vv);
    Jb.multiply(zz, jbzz);
    vv += jbzz;

    // vv_{i+1} = vv_{i+1} - beta * v_{i-1}
    if (jiter > 0) vv.axpy(-beta, *vvecs_[jiter - 1]);

    // alpha_{i} = zz_{i+1}^T vv_{i}
    double alpha = dot_product(zz, vv);

    if (alpha <= 0.0) {
      error_messages_.push_back(
          "SQRTBPLanczos: stopping J'' not positive definite");
      break;
    }

    // vv_{i+1} = vv_{i+1} - alpha_{i} v_{i}
    vv.axpy(-alpha, *vvecs_[jiter]);

    // Re-orthogonalization
    for (size_t jj = 0; jj < jiter; ++jj) {
      double proj = dot_product(vv, *zvecs_[jj]);
      vv.axpy(-proj, *vvecs_[jj]);
    }

    lmp_.inverseMultiply(vv, zz);  // zz = precond vv

    // beta_{i+1} = sqrt( zz_{i+1}^t, vv_{i+1} )
    beta = sqrt(dot_product(zz, vv));

    // v_{i+1} = v_{i+1} / beta_{i+1}
    vv *= 1.0 / beta;
    // zz_{i+1} = zz_{i+1} / beta_{i+1}
    zz *= 1.0 / beta;

    // vvecs[i+1] = v_{i+1}
    vvecs_.push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(vv)));
    ritzPairs_.vVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(vv)));
    // zvecs[i+1] = z_{i+1}
    zvecs_.push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(zz)));
    ritzPairs_.zVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(zz)));

    alphas_.push_back(alpha);
    ritzPairs_.alphas().push_back(alpha);

    if (jiter == 0) {
      ss.push_back(beta0 / alpha);
      dd.push_back(beta0);
    } else {
      // Solve the tridiagonal system T_{i} s_{i} = beta0 * e_1
      dd.push_back(beta0 * dot_product(*zvecs_[0], vv));
      oops::TriDiagSolve(alphas_, betas_, dd, ss);
    }

    betas_.push_back(beta);
    ritzPairs_.betas().push_back(beta);

    // Compute the quadratic cost function
    this->calcQuadCost(rr, ss);  // Check this part !!!!!!!!

    // Compute the Ritz information
    bool ritzFailureCheckFlag = this->calcRitzInformation();

    // Convergence criterion
    if (linfoconv) {
      // Information content based convergence criterion
      //  - Information content defined as 0.5 dv^T dv
      //  - We should measure the change in the solution vector length
      //    in the un-preconditioned metric
      //  - This would require additional call to the preconditioner
      double zdfi_old = zdfi;
      zdfi = costJbLocal_;
      ASSERT(zdfi != 0.0);
      normReduction = (zdfi - zdfi_old) / zdfi;
      oops::Log::info() << "SQRTBPLanczos end of iteration " << jiter + 1 << std::endl
                  << "  Relative information content gain (" << std::setw(2)
                  << jiter + 1 << ") = " << util::full_precision(normReduction)
                  << std::endl;
    } else {
      // Gradient reduction based convergence criterion
      // Gradient norm in precond metric --> sqrt(r'z) --> beta * s_{i}
      double rznorm = beta * std::abs(ss[jiter]);
      normReduction = rznorm / beta0;
      oops::Log::info() << "SQRTBPLanczos end of iteration " << jiter + 1 << std::endl
                        << "  Norm reduction (" << std::setw(2) << jiter + 1
                        << ") = " << util::full_precision(normReduction) << std::endl;
    }

    // Print iteration summary
    this->printQuadCost(jiter);
    this->printRitzInformation(jiter);

    // Check if convergence criteria are met
    if (ritzFailureCheckFlag) break;
    if (normReduction < tolerance && jiter >= miniter) {
      convergenceReached = true;
      break;
    }

    // Increment iteration counter
    ++jiter;
  }

  // Print summary
  oops::Log::info() << "Summary of SQRTBPLanczos solver: " << std::endl
              << " Information based convergence criterion: " << linfoconv
              << std::endl
              << " Number of iterations performed: " << jiter << std::endl
              << " Maximum allowed number of iterations: " << maxiter
              << std::endl
              << " Minimum allowed number of iterations: " << miniter
              << std::endl;
  if (linfoconv) {
    oops::Log::info() << " Required relative gain in information content: "
                << tolerance << std::endl
                << " Last relative gain in information content: "
                << normReduction << std::endl;
  } else {
    oops::Log::info() << " Required reduction in norm of gradient: " << tolerance
                << std::endl
                << " Achieved reduction in norm of gradient: " << normReduction
                << std::endl;
  }
  if (convergenceReached) {
    oops::Log::info() << " Requested convergence criteria met. Well done. "
                << std::endl;
  } else {
    oops::Log::info() << " Failed to meet convergence criteria. " << std::endl;
    for (auto const &err : error_messages_) {
      oops::Log::info() << "  " << err << std::endl;
    }
  }
  oops::Log::info() << std::endl;

  // Calculate the solution
  dv.zero();
  for (size_t jj = 0; jj < ss.size(); ++jj) {
    dv.axpy(ss[jj], *zvecs_[jj]);
  }

  dv += dv0;

  // Process Ritz pairs
  ritzPairs_.process(conf_, "control");

  // Update LMP
  if (oops::SQRTMinimizer<MODEL>::iter_ < oops::SQRTMinimizer<MODEL>::maxiter_)
    lmp_.update(zvecs_, alphas_, betas_);

  // Clean up
  releaseResources();

  return normReduction;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::releaseResources() {
  vvecs_.clear();
  zvecs_.clear();
  alphas_.clear();
  betas_.clear();
  evals_.clear();
  eigvals_.clear();
  evecs_.clear();
  erritz_.clear();
  erritzlm_.clear();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::setupQuadCost(const double &costJ0Jb,
                                                  const double &costJ0JoJc) {
  // J0
  costJ0Jb_ = costJ0Jb;
  costJ0JoJc_ = costJ0JoJc;

  // J0
  costJ0_ = costJ0Jb_ + costJ0JoJc_;

  // reset
  costJ_ = 0;
  costJb_ = 0;
  costJbLocal_ = 0;
  costJoJc_ = 0;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::calcQuadCost(
    const CtrlVec_ &rr, const std::vector<double> &ss) {
  // Compute the quadratic cost function
  // J[du_{i}] = J[0] - 0.5 s_{i}^T Z_{i}^T r_{0}
  // Jb[du_{i}] = 0.5 s_{i}^T V_{i}^T Z_{i} s_{i}
  // Note: here we calculate the quadratic cost function in the preconditioned
  //       space;
  costJ_ = costJ0_;
  costJb_ = costJ0Jb_;
  costJbLocal_ = 0.0;
  for (size_t jj = 0; jj < ss.size(); ++jj) {
    costJ_ -= 0.5 * ss[jj] * dot_product(*zvecs_[jj], rr);
    costJbLocal_ +=
        0.5 * ss[jj] * dot_product(*vvecs_[jj], *zvecs_[jj]) * ss[jj];
  }
  costJb_ += costJbLocal_;
  costJoJc_ = costJ_ - costJb_;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::printQuadCost(const size_t &jiter) const {
  // Print the quadratic cost function
  oops::Log::info() << "  Quadratic cost function: J   (" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJ_) << std::endl
              << "  Quadratic cost function: Jb  (" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJb_) << std::endl
              << "  Quadratic cost function: JoJc(" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJoJc_) << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
bool SQRTBPLanczosEVILMinimizer<MODEL>::calcRitzInformation() {
  ASSERT(alphas_.size() == betas_.size());
  const unsigned nvec = alphas_.size();
  bool ritzFailureCheckFlag = false;
  // Save the leading converged eigenvalue
  std::vector<double> eigvals_old;
  if (eigvals_.size() > 0) eigvals_old.push_back(eigvals_[0]);

  if (nvec > 0) {
    evals_.clear();
    eigvals_.clear();
    evecs_.clear();
    erritz_.clear();
    erritzlm_.clear();
    //  Compute spectrum of tri-diagonal matrix
    oops::TriDiagSpectrum(alphas_, betas_, evals_, evecs_);

    //  Determine the converged eigenvalues
    for (unsigned jiter = 0; jiter < nvec; ++jiter) {
      double erritz = std::abs(evecs_[jiter][nvec - 1] * betas_[nvec - 1]);
      double lambda = evals_[jiter];
      if (lambda < 0) {
        // NEGATIVE RITZ VALUE DETECTED"
        ritzFailureCheckFlag = true;
        error_messages_.push_back("Negative Ritz value detected");
        break;
      }
      //    Add new converged eigenvalues
      double erritzlm = 0.0001 * lambda;
      if (erritz <= erritzlm) {
        eigvals_.push_back(lambda);
        erritz_.push_back(erritz);
        erritzlm_.push_back(erritzlm);
      }
    }

    // Check if Ritz values are not exploding
    if (eigvals_.size() != 0 && eigvals_old.size() != 0) {
      if (eigvals_[0] > 1.01 * eigvals_old[0]) {
        // RITZ VALUES EXPLODE!
        ritzFailureCheckFlag = true;
        oops::Log::info() << "SQRTBPLanczos: Ritz values explode" << std::endl
                    << "Leading Ritz value: " << eigvals_[0] << std::endl
                    << "Leading converged eigenvalue: " << eigvals_old[0]
                    << std::endl;
        error_messages_.push_back("Ritz values explode");
      }
    }
  }
  return ritzFailureCheckFlag;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::printRitzInformation(
    const size_t &jiter) const {
  ASSERT(alphas_.size() == betas_.size());
  const unsigned nvec = alphas_.size();

  if (nvec > 0) {
    oops::Log::info() << "  Converged Ritz values (" << jiter + 1
                << "):" << std::endl;
    for (auto const &eigval : eigvals_) {
      oops::Log::info() << "    Ritz value: " << util::full_precision(eigval)
                  << std::endl;
    }
    for (auto const &err : erritz_) {
      oops::Log::info() << "    Error bounds: " << util::full_precision(err)
                  << std::endl;
    }
    for (auto const &errlm : erritzlm_) {
      oops::Log::info() << "    Error bound limits: " << util::full_precision(errlm)
                  << std::endl;
    }
    oops::Log::info() << std::endl;
  }
  return;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::checkpointLMP(
    eckit::LocalConfiguration &conf) const {
  ASSERT(conf.has("preconditioner"));

  // Checkpoint lmp
  lmp_.checkpoint(conf);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTBPLanczosEVILMinimizer<MODEL>::restartLMP(
    const eckit::Configuration &conf) {
  ASSERT(conf.has("preconditioner"));

  // Restart lmp
  lmp_.restart(conf);
}

// -----------------------------------------------------------------------------

}  // namespace saber
