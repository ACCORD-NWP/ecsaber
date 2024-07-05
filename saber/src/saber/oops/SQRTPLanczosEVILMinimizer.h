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

#include "util/LogbookWriter.h"

#include "util/dot_product.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SQRTPLanczosEVILMinimizer : public oops::SQRTMinimizer<MODEL> {
  using CostFct_ = oops::CostFunction<MODEL>;
  using CtrlVec_ = oops::ControlVector<MODEL>;
  using Jbmat_ = oops::JbMatrix<MODEL>;
  using UtHtRinvHU_ = oops::UtHtRinvHUMatrix<MODEL>;

 public:
  const std::string classname() const final
    { return "SQRTPLanczosEVILMinimizer"; }

  /// Constructor, destructor
  SQRTPLanczosEVILMinimizer(const eckit::Configuration &,
                        const CostFct_ &);
  ~SQRTPLanczosEVILMinimizer() {}

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
  void calcQuadCost(const CtrlVec_ &, const CtrlVec_ &,
                    const std::vector<double> &);
  void printQuadCost(const size_t &) const;

  /// Other diagnostics
  bool calcRitzInformation();
  void printRitzInformation(const size_t &) const;

  /// Data memebers
  const eckit::LocalConfiguration conf_;
  const CostFct_ &J_;

  /// LMP
  oops::SpectralSqrtLMP<MODEL> lmp_;

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
  double costJbLocal_;  // component of costJb from current minimization
  double costJoJc_;
};

// =============================================================================

template <typename MODEL>
SQRTPLanczosEVILMinimizer<MODEL>::SQRTPLanczosEVILMinimizer(
    const eckit::Configuration &conf, const CostFct_ &J)
    : oops::SQRTMinimizer<MODEL>(conf, J),
      conf_(conf),
      J_(J),
      lmp_(conf_, J),
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
      costJbLocal_(0),
      costJoJc_(0) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
double SQRTPLanczosEVILMinimizer<MODEL>::solve(
    CtrlVec_ &dv, CtrlVec_ &rr, const Jbmat_ &Jb, const UtHtRinvHU_ &UtHtRinvHU,
    const bool &linfoconv, const size_t &maxiter, const size_t &miniter,
    const double &tolerance, const double &costJ0Jb, const double &costJ0JoJc) {
  // dv   control vector

  std::vector<double> ss;
  std::vector<double> dd;

  // J0
  this->setupQuadCost(costJ0Jb, costJ0JoJc);
  oops::Log::info() << "Init Quadratic cost J0 = " << costJ0_ << std::endl;

  // Change resolution of LMP vectors
  lmp_.changeResolution();

  // Postprocess preconditioning vectors
  lmp_.postprocess();

  // v_{0} = P^-1/2 r_{0}
  std::unique_ptr<CtrlVec_> mrr(new CtrlVec_(rr, false));
  lmp_.inverseMultiplySqrt(rr, *mrr);
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
  while (jiter < maxiter) {
    size_t jiterp1 = jiter + 1;
    oops::Log::info() << "SQRTPLanczosEVIL Starting Iteration " << jiterp1 << std::endl;

    util::LogbookWriter<size_t> log("InnerLoop", jiterp1);

    // v_{i+1} = P^-1/2 ( I + U^T H^T R^{-1} H U ) P^-1/2 v_{i}
    lmp_.inverseMultiplySqrt(*vv, *mvv);

    UtHtRinvHU.multiply(*mvv, *zz);

    Jb.multiply(*mvv, *jbzz);

    *zz += *jbzz;

    lmp_.inverseMultiplySqrt(*zz, *vv);

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
      oops::TriDiagSolve(ritzPairs_.alphas(), ritzPairs_.betas(), dd, ss);
    }

    ritzPairs_.betas().push_back(beta);

    // Compute the quadratic cost function
    this->calcQuadCost(dv, *mrr, ss);

    // Compute the Ritz information
    bool ritzErrorDetected = this->calcRitzInformation();

    // Gradient norm in precond metric --> sqrt(r'z) --> beta * s_{i}
    double gradNorm = beta * std::abs(ss[jiter]);
    double gradNormReduction = gradNorm / beta0;

    oops::Log::info() << "SQRTPLanczosEVIL end of iteration " << jiterp1 << std::endl
                      << "  Gradient norm (" << std::setw(2) << jiterp1
                      << ") = " << util::full_precision(gradNorm) << std::endl
                      << "  Gradient norm reduction (" << std::setw(2) << jiterp1
                      << ") = " << util::full_precision(gradNormReduction)
                      << std::endl;

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
      oops::Log::info() << "  Relative information content gain (" << std::setw(2)
                        << jiterp1 << ") = " << util::full_precision(normReduction)
                        << std::endl;
    } else {
      // Gradient reduction based convergence criterion
      normReduction = gradNormReduction;
    }
    if (normReduction < tolerance && jiter >= miniter)
      convergenceReached = true;

    util::LogbookWriter<double> log_norm("NormReduction", normReduction);

    // Print iteration summary
    this->printQuadCost(jiter);
    this->printRitzInformation(jiter);

    // Increment iteration counter
    ++jiter;

    // Compute online diagnostics
    UtHtRinvHU.computeDiagnostics(ss, convergenceReached || (jiter == maxiter));

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
  oops::Log::info() << "Summary of SQRTPLanczosEVIL solver: " << std::endl
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
    //  Note that soft errors do not result in the application being aborted;
    // minimization terminates once a soft error is detected, the increment
    // is still computed.
    oops::Log::info() << " Failed to meet convergence criteria. " << std::endl;
    for (auto const &err : soft_error_messages_) {
      oops::Log::info() << "  " << err << std::endl;
    }
  }
  oops::Log::info() << std::endl;

  // Calculate the solution
  std::unique_ptr<CtrlVec_> dv0(new CtrlVec_(dv));
  std::unique_ptr<CtrlVec_> mdv(new CtrlVec_(dv, false));
  mdv->zero();
  for (size_t jj = 0; jj < ss.size(); ++jj) {
    mdv->axpy(ss[jj], ritzPairs_.vVEC(jj));
  }

  // Transform control variable back to space with scalar product
  lmp_.inverseMultiplySqrt(*mdv, dv);

  dv += *dv0;

  // Free memory
  dv0.reset();
  mdv.reset();

  // Process Ritz pairs
  ritzPairs_.process(conf_);

  // Update LMP
  if (oops::SQRTMinimizer<MODEL>::iter_ < oops::SQRTMinimizer<MODEL>::maxiter_)
    lmp_.update(ritzPairs_.vVEC(), ritzPairs_.alphas(), ritzPairs_.betas());

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
  costJbLocal_ = 0;
  costJoJc_ = 0;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::calcQuadCost(const CtrlVec_ &dv0,
                                                const CtrlVec_ &rr,
                                                const std::vector<double> &ss) {
  // Compute the quadratic cost function
  // J[du_{i}] = J[0] - 0.5 s_{i}^T Z_{i}^T r_{0}
  // Jb[du_{i}] = 0.5 s_{i}^T V_{i}^T Z_{i} s_{i}
  // Note: here we calculate the quadratic cost function in the preconditioned
  //       space;
  costJ_ = costJ0_;
  costJb_ = costJ0Jb_;
  costJbLocal_ = 0.0;
  for (size_t jj = 0; jj < ss.size(); ++jj) {
    costJ_ -= 0.5 * ss[jj] * dot_product(ritzPairs_.vVEC(jj), rr);
    costJbLocal_ += 0.5 * ss[jj] *
                    dot_product(ritzPairs_.vVEC(jj), ritzPairs_.vVEC(jj)) *
                    ss[jj];
    /* EA TODO:
     * if (not lclassicincrform) {
     *  costJbLocal_ += ss[jj] * dot_product(ritzPairs_.vVEC(jj), dv0);
    }*/
  }
  costJb_ += costJbLocal_;
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
  oops::Log::info() << "  Quadratic cost function: J   (" << std::setw(2) << jiter + 1
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
  const unsigned nvec = ritzPairs_.alphas().size();
  bool ritzErrorDetected = false;

  if (nvec > 0) {
    eig_.clear();
    erreig_.clear();
    erreiglm_.clear();
    std::vector<double> ritzvals;
    std::vector<std::vector<double> > ritzvecs;
    //  Compute spectrum of tri-diagonal matrix
    oops::TriDiagSpectrum(ritzPairs_.alphas(), ritzPairs_.betas(), ritzvals,
                          ritzvecs);

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
      double erritz =
          std::abs(ritzvecs[jiter][nvec - 1] * ritzPairs_.betas()[nvec - 1]);
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
        oops::Log::info() << "SQRTPLanczosEVIL: Ritz values explode" << std::endl
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
  return ritzErrorDetected;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::printRitzInformation(
    const size_t &jiter) const {
  ASSERT(ritzPairs_.alphas().size() == ritzPairs_.betas().size());
  const unsigned nvec = ritzPairs_.alphas().size();

  if (nvec > 0) {
    oops::Log::info() << "  Converged Ritz values (" << jiter + 1
                      << "):" << std::endl;
    for (auto const &eigval : eig_) {
      oops::Log::info() << "    Ritz value: " << util::full_precision(eigval)
                        << std::endl;
    }
    for (auto const &err : erreig_) {
      oops::Log::info() << "    Error bounds: " << util::full_precision(err)
                        << std::endl;
    }
    for (auto const &errlm : erreiglm_) {
      oops::Log::info() << "    Error bound limits: " << util::full_precision(errlm)
                        << std::endl;
    }
    oops::Log::info() << std::endl;
  }
  return;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::checkpointLMP(
    eckit::LocalConfiguration &conf) const {
  ASSERT(conf.has("preconditioner"));

  // Checkpoint lmp
  lmp_.checkpoint(conf);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SQRTPLanczosEVILMinimizer<MODEL>::restartLMP(
    const eckit::Configuration &conf) {
  ASSERT(conf.has("preconditioner"));

  // Restart lmp
  lmp_.restart(conf);
}

// -----------------------------------------------------------------------------

}  // namespace saber
