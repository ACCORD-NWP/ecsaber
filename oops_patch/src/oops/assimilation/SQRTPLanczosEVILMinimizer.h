/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/assimilation/ControlVector.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/JbMatrix.h"
#include "oops/assimilation/RitzPairs.h"
#include "oops/assimilation/SQRTMinimizer.h"
#include "oops/assimilation/SpectralSqrtLMP.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/assimilation/TriDiagSpectrum.h"
#include "oops/assimilation/UtHtRinvHUMatrix.h"
#include "util/LogbookWriter.h"
#include "util/Logger.h"
#include "util/abor1_cpp.h"
#include "util/dot_product.h"
#include "util/formats.h"

namespace oops {

/// SQRTPLanczosMinimizer
/*!
 * \brief Preconditioned Lanczos solver for the square root preconditioned
 * formulation of the Hessian
 *
 * The Hessian must be square, symmetric and positive definite.

 * On entry:
 * -    A       = \f$ (U^T B^{-1} U + U^T H^T R^{-1} U) \f$
 * -    dv      =  starting point, \f$ dv_0 \f$.
 * -    rr      = \f$ (-\sum dv^{b}_{i} + ) U^T H^T R^{-1} d  - A dv_0 \f$
 * - or
 * -    rr      = \f$ U^T H^T R^{-1} [d + H (x_b - x_{fg})] - A dv_0 \f$

 * On exit, dv will contain the solution such that \f$ dx = U dv \f$ or
 * \f$ dx = U dv + x_{b} - x_{fg} \f$
 *  The solution is recovered in the SQRTMinimizer class
 *  The return value is the achieved reduction in preconditioned residual
 *  norm/information content gain.
 *
 *  Iteration will stop if the maximum iteration limit "maxIter" is reached
 *  or if the residual norm reduces by a factor of "tolerance"/information
 *  content based convergence criterion is reached.
 *
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class SQRTPLanczosEVILMinimizer : public SQRTMinimizer<MODEL> {
  using CostFct_ = CostFunction<MODEL>;
  using CtrlVec_ = ControlVector<MODEL>;
  using Jbmat_ = JbMatrix<MODEL>;
  using UtHtRinvHU_ = UtHtRinvHUMatrix<MODEL>;

 public:
  const std::string classname() const final { return "SQRTPLanczosEVILMinimizer"; }

  /// Constructor, destructor
  SQRTPLanczosEVILMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~SQRTPLanczosEVILMinimizer() {}

 private:
  double solve(CtrlVec_ &, CtrlVec_ &, const Jbmat_ &,
               const UtHtRinvHU_ &) final;

  void releaseResources();

  /// Checkpointing, restart of the limited memory preconditioner
  void checkpointLMP(eckit::LocalConfiguration &) const final;
  void restartLMP(const eckit::Configuration &) final;

  /// Quadratic cost function calculations
  void setupQuadCost(const double &costJ0Jb, const double &costJ0JoJc);
  void calcQuadCost(const CtrlVec_ &, const CtrlVec_ &, const CtrlVec_ &,
                    const CtrlVec_ &, const CtrlVec_ &, const bool &);
  void printQuadCost(const size_t &) const;

  /// Other diagnostics
  bool calcRitzInformation();
  void printRitzInformation(const size_t &) const;

  /// Data memebers
  const eckit::LocalConfiguration conf_;
  const CostFct_ &J_;

  /// LMP
  std::unique_ptr<SpectralSqrtLMP<MODEL>> lmp_;

  /// Local
  size_t itheta1_{0};
  double ztheta1_{0.0};
  RitzPairs<CtrlVec_> ritzPairs_;

  std::vector<double> eig_;
  std::vector<double> erreig_;
  std::vector<double> erreiglm_;

  std::vector<std::string> soft_error_messages_;

  /// Quadratic cost function calculations data members
  double costJ0_;
  double costJ0Jb_;
  double costJ0JoJc_;

  double costJ_;
  double costJm1_;
  double costJb_;
  double costJbCurrentMin_;  // component of costJb from current minimization
  double costJbCurrentMinPrecSpace_;  // component of costJb from current
                                      // minimization in preconditioned space
  double costJoJc_;
};

// =============================================================================

template <typename MODEL>
SQRTPLanczosEVILMinimizer<MODEL>::SQRTPLanczosEVILMinimizer(
    const eckit::Configuration &conf, const CostFct_ &J)
    : SQRTMinimizer<MODEL>(conf, J),
      conf_(conf),
      J_(J),
      lmp_(),
      eig_(),
      erreig_(),
      erreiglm_(),
      soft_error_messages_(),
      costJ0_(0),
      costJ0Jb_(0),
      costJ0JoJc_(0),
      costJ_(0),
      costJm1_(0),
      costJb_(0),
      costJbCurrentMin_(0),
      costJbCurrentMinPrecSpace_(0),
      costJoJc_(0) {
  Log::trace() << classname() << "::SQRTPLanczosMinimizer() starting"
               << std::endl;
  eckit::LocalConfiguration precondConf;
  if (conf.has("preconditioner")) {
    precondConf = conf.getSubConfiguration("preconditioner");
  }
  lmp_ = std::move(SqrtLMPFactory<MODEL>::create(precondConf, J));
  Log::trace() << classname() << "::SQRTPLanczosMinimizer() done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double SQRTPLanczosEVILMinimizer<MODEL>::solve(CtrlVec_ &dv, CtrlVec_ &rr,
                                               const Jbmat_ &Jb,
                                               const UtHtRinvHU_ &UtHtRinvHU) {
  // Setup controls
  const bool &linfoConv = SQRTMinimizer<MODEL>::linfoConv_;
  const bool &lclassicChaVar = SQRTMinimizer<MODEL>::lclassicChaVar_;
  const size_t &maxIter = SQRTMinimizer<MODEL>::ninner_;
  const size_t &minIter = SQRTMinimizer<MODEL>::ninnerMin_;
  const double &tolerance = SQRTMinimizer<MODEL>::gnreduc_;
  const double &costJ0Jb = SQRTMinimizer<MODEL>::costJ0Jb_;
  const double &costJ0JoJc = SQRTMinimizer<MODEL>::costJ0JoJc_;
  const CtrlVec_ &gradJb = *(SQRTMinimizer<MODEL>::gradJb_);
  const size_t &outerIter = SQRTMinimizer<MODEL>::outerIteration_;

  // Auxiliary vectors
  std::vector<double> ss;
  std::vector<double> dd;

  // Initial quadcost = J0
  this->setupQuadCost(costJ0Jb, costJ0JoJc);
  Log::info() << "Init Quadratic cost J0 = " << costJ0_ << std::endl;

  // Change resolution of LMP vectors
  lmp_->changeResolution();

  // Postprocess preconditioning vectors.
  // If the LMP is based on randomization, this will generate LMP vectors.
  // The initial residual is rr only if dv = 0
  lmp_->postprocess(UtHtRinvHU, Jb, rr, outerIter);

  // v_{0} = P^-1/2 r_{0}
  std::unique_ptr<CtrlVec_> mrr(new CtrlVec_(rr, false));
  lmp_->inverseMultiplySqrt(rr, *mrr);

  std::unique_ptr<CtrlVec_> vv(new CtrlVec_(*mrr));

  // beta_{0} = sqrt( v_{0}^T v_{0} )
  double beta = sqrt(dot_product(*vv, *vv));
  const double beta0 = beta;

  // v_{1} = r_{0} / beta_{0}
  *vv *= 1.0 / beta;

  // vvecs[0] = v_{1} ---> for re-orthogonalization
  ritzPairs_.vVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(*vv)));

  double normReduction = 1.0;
  double zdfi = 0.0;
  bool convergenceReached = false;
  size_t jiter = 0;

  std::unique_ptr<CtrlVec_> mvv(new CtrlVec_(dv, false));
  std::unique_ptr<CtrlVec_> jbzz(new CtrlVec_(dv, false));
  std::unique_ptr<CtrlVec_> zz(new CtrlVec_(rr, false));
  std::unique_ptr<CtrlVec_> mds(new CtrlVec_(dv, false));
  std::unique_ptr<CtrlVec_> ds(new CtrlVec_(dv, false));
  while (jiter < maxIter) {
    size_t jiterp1 = jiter + 1;
    Log::info() << "SQRTPLanczosEVIL Starting Iteration " << jiterp1 << std::endl;

    util::LogbookWriter<size_t> log("InnerLoop", jiterp1);

    // v_{i+1} = P^-1/2 ( I + U^T H^T R^{-1} H U ) P^-1/2 v_{i}
    lmp_->inverseMultiplySqrt(*vv, *mvv);

    UtHtRinvHU.multiply(*mvv, *zz);

    Jb.multiply(*mvv, *jbzz);

    *zz += *jbzz;

    lmp_->inverseMultiplySqrt(*zz, *vv);

    // v_{i+1} = v_{i+1} - beta * v_{i-1}
    if (jiter > 0) vv->axpy(-beta, ritzPairs_.vVEC(jiter - 1));

    // alpha_{i} = v_{i+1}^T v_{i}
    double alpha = dot_product(ritzPairs_.vVEC(jiter), *vv);
    if (alpha <= 0.0) {
      soft_error_messages_.push_back(
          "SQRTPLanczosEVIL: stopping J'' not positive definite");
      break;
    }

    // v_{i+1} = v_{i+1} - alpha_{i} v_{i}
    vv->axpy(-alpha, ritzPairs_.vVEC(jiter));

    // Re-orthogonalization
    for (size_t jj = 0; jj < jiter; ++jj) {
      double proj = dot_product(*vv, ritzPairs_.vVEC(jj));
      vv->axpy(-proj, ritzPairs_.vVEC(jj));
    }

    // beta_{i+1} = sqrt( vv_{i+1}^t, vv_{i+1} )
    beta = sqrt(dot_product(*vv, *vv));

    // v_{i+1} = v_{i+1} / beta_{i+1}
    *vv *= 1.0 / beta;

    // vvecs[i+1] = v_{i+1}
    ritzPairs_.vVEC().push_back(std::shared_ptr<CtrlVec_>(new CtrlVec_(*vv)));

    ritzPairs_.alphas().push_back(alpha);

    if (jiter == 0) {
      ss.push_back(beta0 / alpha);
      dd.push_back(beta0);
    } else {
      // Solve the tridiagonal system T_{i} s_{i} = beta0 * e_1
      dd.push_back(beta0 * dot_product(ritzPairs_.vVEC(0), *vv));
      TriDiagSolve(ritzPairs_.alphas(), ritzPairs_.betas(), dd, ss);
    }

    ritzPairs_.betas().push_back(beta);

    // Reconstruct the control variable in the B preconditioned space
    mds->zero();
    for (size_t jj = 0; jj < ss.size(); ++jj) {
      mds->axpy(ss[jj], ritzPairs_.vVEC(jj));
    }

    // Transform the control variable back to B preconditioned space
    lmp_->inverseMultiplySqrt(*mds, *ds);

    // Compute the quadratic cost function
    this->calcQuadCost(dv, *mds, *ds, rr, gradJb, lclassicChaVar);

    // Compute the Ritz information
    const bool ritzErrorDetected = this->calcRitzInformation();

    // Gradient norm in precond metric --> sqrt(r'z) --> beta * s_{i}
    double gradNorm = beta * std::abs(ss[jiter]);
    double gradNormReduction = gradNorm / beta0;

    Log::info() << "SQRTPLanczosEVIL end of iteration " << jiterp1 << std::endl
                << "  Gradient norm (" << std::setw(2) << jiterp1
                << ") = " << util::full_precision(gradNorm) << std::endl
                << "  Gradient norm reduction (" << std::setw(2) << jiterp1
                << ") = " << util::full_precision(gradNormReduction)
                << std::endl;

    // Convergence criterion
    if (linfoConv) {
      // Information content based convergence criterion
      //  - Information content defined as 0.5 dv^T dv
      //  - We should measure the change in the solution vector length
      //    in the un-preconditioned metric
      //  - This would require additional call to the preconditioner
      double zdfi_old = zdfi;
      zdfi = costJbCurrentMinPrecSpace_;
      ASSERT(zdfi != 0.0);
      normReduction = (zdfi - zdfi_old) / zdfi;
      Log::info() << "  Relative information content gain (" << std::setw(2)
                  << jiterp1 << ") = " << util::full_precision(normReduction)
                  << std::endl;
    } else {
      // Gradient reduction based convergence criterion
      normReduction = gradNormReduction;
    }
    if (normReduction < tolerance && jiter >= minIter)
      convergenceReached = true;

    util::LogbookWriter<double> log_norm("NormReduction", normReduction);

    // Print iteration summary
    this->printQuadCost(jiter);
    this->printRitzInformation(jiter);

    // Increment iteration counter
    ++jiter;

    // Compute online diagnostics
    UtHtRinvHU.computeDiagnostics(ss, convergenceReached || (jiter == maxIter));

    // Check if convergence criteria are met
    if (convergenceReached || ritzErrorDetected) break;
  }
  // Free memory
  mrr.reset();
  vv.reset();
  mvv.reset();
  jbzz.reset();
  zz.reset();

  // Print summary
  Log::info() << "Summary of SQRTPLanczos solver: " << std::endl
              << " Information based convergence criterion: " << linfoConv
              << std::endl
              << " Number of iterations performed: " << jiter << std::endl
              << " Maximum allowed number of iterations: " << maxIter
              << std::endl
              << " Minimum allowed number of iterations: " << minIter
              << std::endl;
  if (linfoConv) {
    Log::info() << " Required relative gain in information content: "
                << tolerance << std::endl
                << " Last relative gain in information content: "
                << normReduction << std::endl;
  } else {
    Log::info() << " Required reduction in norm of gradient: " << tolerance
                << std::endl
                << " Achieved reduction in norm of gradient: " << normReduction
                << std::endl;
  }
  if (convergenceReached) {
    Log::info() << " Requested convergence criteria met. Well done. "
                << std::endl;
  } else {
    //  Note that soft errors do not result in the application being aborted;
    // minimization terminates once a soft error is detected, the increment
    // is still computed.
    Log::info() << " Failed to meet convergence criteria. " << std::endl;
    for (auto const &err : soft_error_messages_) {
      Log::info() << "  " << err << std::endl;
    }
  }
  Log::info() << std::endl;

  // Update the solution
  dv += *ds;

  // Free memory
  mds.reset();

  // Process Ritz pairs
  ritzPairs_.process(conf_, "control");

  // Update LMP
  if (SQRTMinimizer<MODEL>::outerIteration_ <
      SQRTMinimizer<MODEL>::lastOuterIteration_)
    lmp_->update(ritzPairs_.vVEC(), ritzPairs_.alphas(), ritzPairs_.betas());

  // Clean up
  releaseResources();

  return normReduction;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::releaseResources() {
  ritzPairs_.vVEC().clear();
  ritzPairs_.alphas().clear();
  ritzPairs_.betas().clear();
  eig_.clear();
  erreig_.clear();
  erreiglm_.clear();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::setupQuadCost(const double &costJ0Jb,
                                                     const double &costJ0JoJc) {
  // J0
  costJ0Jb_ = costJ0Jb;
  costJ0JoJc_ = costJ0JoJc;

  // J0
  costJ0_ = costJ0Jb_ + costJ0JoJc_;

  // reset
  costJ_ = 0;
  costJm1_ = 0;
  costJb_ = 0;
  costJbCurrentMin_ = 0;
  costJbCurrentMinPrecSpace_ = 0;
  costJoJc_ = 0;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::calcQuadCost(
    const CtrlVec_ &dv0, const CtrlVec_ &mds, const CtrlVec_ &ds,
    const CtrlVec_ &rr, const CtrlVec_ &gradJb, const bool &lclassicChaVar) {
  // Compute the quadratic cost function in the state space
  //
  // At ith iteration of SQRTPLanczos, the solution is
  // dv_{i} = dv_{0}  + P V_{i} s_{i} where P is a second-level preconditioner
  //
  // We can calculate the quadratic cost function as:
  // J[dv_{i}] = J[dv_{0}] - 0.5 r_{0}^T P V_{i} s_{i} with
  //     r_{0} = U^T b - U^T A U dv_{0}
  //
  // If lclassicChaVar = true,
  //    dx = U dv
  //    dv_{0} = 0
  //    Jb[dv_{i}] = Jb[dv_{0}] + 0.5 dv_{i}^T U^T Binv U dv_{i} - gradJb^T
  //    dv_{i}
  // else
  //    dx = U dv - xk + xb
  //    dv_{0}^{k} = dv_{j}^{k} with j being the last iteration of the previous
  //    (k)-th system. Warm start initial point. Jb[dv_{i}] = dv_{i}^T U^T Binv
  //    B dv_{i}
  //
  // Note that for two cases J[dv_{0}] are the same which is equivalent to
  // J[dx = 0]. Note that U^T Binv B is assumed to be an identity matrix.

  // Initialize cost function values J[dv_{0}] and Jb[dv_{0}]
  costJ_ = costJ0_;
  costJb_ = costJ0Jb_;

  // Calculate the quadratic cost function J[dv_{i}]
  costJ_ -= 0.5 * dot_product(rr, ds);

  // Calculate Jb part of the quadratic cost function: Jb[dv_{i}]
  if (lclassicChaVar) {
    // dv = 0 + ds = ds (zero initial guess)
    // SG: This diagnostics should be removed for later versions.
    costJbCurrentMinPrecSpace_ =
        dot_product(mds, gradJb) + 0.5 * dot_product(mds, mds);
    costJbCurrentMin_ = dot_product(ds, gradJb) + 0.5 * dot_product(ds, ds);
    costJb_ += costJbCurrentMin_;
  } else {
    // dv = dv_0 + ds
    // MC : switch to relying on costJbCurrentMin beyond CY49R1
    costJbCurrentMinPrecSpace_ = 0.5 * dot_product(mds, mds);
    costJbCurrentMin_ = 0.5 * dot_product(ds, ds);
    costJb_ = costJbCurrentMin_ +
              0.5 * (dot_product(dv0, dv0) + 2 * dot_product(ds, dv0));
  }

  // Calculate Jo part of the quadratic cost function:
  // // Jo[dv_{i}] + Jc[dv_{i}] = J[dv_{i}] - Jb[dv_{i}]
  costJoJc_ = costJ_ - costJb_;

  // Check for divergence

  // 1. Check if the quadratic cost function is positive
  if (costJ_ < 0)
    ABORT("SQRTPLanczosEVIL: Fatal error detected: negative quadratic J");

  // 2. Check if the quadratic cost function is decreasing monotonically
  if (costJm1_ > 0 && costJ_ > costJm1_)
    ABORT("SQRTPLanczosEVIL: Fatal error detected: growing quadratic J");

  costJm1_ = costJ_;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::printQuadCost(const size_t &jiter) const {
  // Print the quadratic cost function
  Log::info() << "  Quadratic cost function: J   (" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJ_) << std::endl
              << "  Quadratic cost function: Jb  (" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJb_) << std::endl
              << "  Quadratic cost function: JoJc(" << std::setw(2) << jiter + 1
              << ") = " << util::full_precision(costJoJc_) << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
bool SQRTPLanczosEVILMinimizer<MODEL>::calcRitzInformation() {
  ASSERT(ritzPairs_.alphas().size() == ritzPairs_.betas().size());
  const int nvec = ritzPairs_.alphas().size();
  bool ritzErrorDetected = false;

  if (nvec > 0) {
    eig_.clear();
    erreig_.clear();
    erreiglm_.clear();
    std::vector<double> ritzvals;
    std::vector<std::vector<double>> ritzvecs;
    //  Compute spectrum of tri-diagonal matrix
    oops::TriDiagSpectrum(ritzPairs_.alphas(), ritzPairs_.betas(), ritzvals, ritzvecs);

    std::vector<double> erritzvals;
    std::vector<double> erritzlmvals;
    //  Compute error bounds and order Ritz values (largest to smallest)
    for (int jiter = nvec - 1; jiter >= 0; --jiter) {
      double lambda = ritzvals[jiter];
      if (lambda < 0) {
        //  NEGATIVE RITZ VALUE DETECTED"
        ritzErrorDetected = true;
        soft_error_messages_.push_back("Negative Ritz value detected");
        break;
      }
      double erritz = std::abs(ritzvecs[jiter][nvec - 1] * ritzPairs_.betas()[nvec - 1]);
      double erritzlm = 0.0001 * lambda;
      //  Store Ritz values and their error bounds
      erritzvals.push_back(erritz);
      erritzlmvals.push_back(erritzlm);
    }
    std::reverse(ritzvals.begin(), ritzvals.end());

    // Leading eigenvalue explosion test
    if (eig_.size() > 0) {
      if (ritzvals[itheta1_] > 1.01 * ztheta1_) {
        // RITZ VALUES EXPLODE!
        ritzErrorDetected = true;
        Log::info() << "SQRTPLanczosEVIL: Ritz values explode" << std::endl
                    << "Leading Ritz value: " << ritzvals[itheta1_] << std::endl
                    << "Leading converged eigenvalue: " << ztheta1_
                    << std::endl;
        soft_error_messages_.push_back("Ritz values explode");
      }
    }

    for (int jiter = 0; jiter <= nvec - 1; ++jiter) {
      //  Store converged eigenvalues and their error bounds (largest to
      //  smallest)
      if (erritzvals[jiter] < erritzlmvals[jiter]) {
        eig_.push_back(ritzvals[jiter]);
        erreig_.push_back(erritzvals[jiter]);
        erreiglm_.push_back(erritzlmvals[jiter]);
      }
    }

    for (int jiter = nvec - 1; jiter >= 0; --jiter) {
      //  Save leading converged eigenvalue
      if (erritzvals[jiter] <= erritzlmvals[jiter]) {
        ztheta1_ = ritzvals[jiter];
        itheta1_ = jiter;
      }
    }
  }
  Log::trace() << classname() << "::calcRitzInformation() done" << std::endl;
  return ritzErrorDetected;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::printRitzInformation(
    const size_t &jiter) const {
  ASSERT(ritzPairs_.alphas().size() == ritzPairs_.betas().size());
  const unsigned nvec = ritzPairs_.alphas().size();

  if (nvec > 0) {
    Log::info() << "  Converged Ritz values (" << jiter + 1
                << "):" << std::endl;
    for (auto const &eigval : eig_) {
      Log::info() << "    Ritz value: " << util::full_precision(eigval)
                  << std::endl;
    }
    for (auto const &err : erreig_) {
      Log::info() << "    Error bounds: " << util::full_precision(err)
                  << std::endl;
    }
    for (auto const &errlm : erreiglm_) {
      Log::info() << "    Error bound limits: " << util::full_precision(errlm)
                  << std::endl;
    }
    Log::info() << std::endl;
  }
  return;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::checkpointLMP(
    eckit::LocalConfiguration &conf) const {
  lmp_->checkpoint(conf);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::restartLMP(
    const eckit::Configuration &conf) {
  lmp_->restart(conf);
}

// -----------------------------------------------------------------------------

}  // namespace oops
