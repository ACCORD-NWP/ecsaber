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
#include <memory>
#include <string>
#include <vector>

#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DRMinimizer.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/assimilation/RitzPairs.h"
#include "oops/assimilation/SpectralLMP.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "util/Logger.h"
#include "util/dot_product.h"
#include "util/formats.h"

namespace oops {

/// DRPLanczos Minimizer
/*!
 * \brief Derber-Rosati Preconditioned Lanczos solver.
 *
 * This solver is the Lanczos version of the DRPCG algorithm
 * It solves \f$ Ax=b\f$ for the particular case \f$ A=B^{-1}+C\f$,
 * without requiring the application of \f$ B^{-1}\f$.
 *
 * A must be square, symmetric, positive definite.
 *
 * A preconditioner must be supplied that, given a vector q, returns an
 * approximation to \f$ (AB)^{-1} q\f$. Possible preconditioning
 * is detailed in S. Gurol, PhD Manuscript, 2013.
 * Note that the traditional \f$ B\f$-preconditioning corresponds to
 * precond=\f$I\f$.
 *
 * On entry:
 * -    dx      = starting point.
 * -    dxh     = starting point, \f$ B^{-1} dx_{0}\f$.
 * -    rr      = residual at starting point.
 * -    B       = \f$ B \f$.
 * -    C       = \f$ C \f$.
 * -    precond = preconditioner \f$ F_k \approx (AB)^{-1} \f$.
 *
 * On exit, dxh will contain \f$ B^{-1} x\f$ where x is the solution.
 * The return value is the achieved reduction in residual norm.
 *
 * Iteration will stop if the maximum iteration limit "maxiter" is reached
 * or if the residual norm reduces by a factor of "tolerance".
 *
 * Each matrix must implement a method:
 * - void multiply(const VECTOR&, VECTOR&) const
 *
 * which applies the matrix to the first argument, and returns the
 * matrix-vector product in the second. (Note: the const is optional, but
 * recommended.)
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class DRPLanczosEVILMinimizer : public DRMinimizer<MODEL> {
  using Bmat_ = BMatrix<MODEL>;
  using CostFct_ = CostFunction<MODEL>;
  using CtrlInc_ = ControlIncrement<MODEL>;
  using HtRinvH_ = HtRinvHMatrix<MODEL>;

 public:
  const std::string classname() const override { return "DRPLanczosEVILMinimizer"; }
  DRPLanczosEVILMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~DRPLanczosEVILMinimizer() {}

 private:
  double solve(CtrlInc_ &, CtrlInc_ &, CtrlInc_ &, const Bmat_ &,
               const HtRinvH_ &, const CtrlInc_ &, const double, const double,
               const int, const double) override;

  const eckit::LocalConfiguration config_;
  SpectralLMP<CtrlInc_> lmp_;
};

// =============================================================================

template <typename MODEL>
DRPLanczosEVILMinimizer<MODEL>::DRPLanczosEVILMinimizer(
    const eckit::Configuration &conf, const CostFct_ &J)
    : DRMinimizer<MODEL>(J),
      config_(conf),
      lmp_(conf) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
double DRPLanczosEVILMinimizer<MODEL>::solve(
    CtrlInc_ &dx, CtrlInc_ &dxh, CtrlInc_ &rr, const Bmat_ &B,
    const HtRinvH_ &HtRinvH, const CtrlInc_ &gradJb, const double costJ0Jb,
    const double costJ0JoJc, const int maxiter, const double tolerance) {
  // dx   increment
  // dxh  B^{-1} dx
  // rr   (sum B^{-1} dx_i^{b} +) G^T H^{-1} d

  RitzPairs<CtrlInc_> ritzPairs;
  CtrlInc_ zz(dxh);
  CtrlInc_ pr(dxh);
  CtrlInc_ vv(rr);
  CtrlInc_ mm(rr);

  std::vector<double> ss;
  std::vector<double> dd;

  // J0
  const double costJ0 = costJ0Jb + costJ0JoJc;

  lmp_.update(ritzPairs.vVEC(), ritzPairs.zVEC(), ritzPairs.tVEC(), ritzPairs.alphas(),
    ritzPairs.betas());

  // z_{0} = B LMP r_{0}
  lmp_.multiply(vv, pr);
  B.multiplyB(pr, zz);

  // beta_{0} = sqrt( z_{0}^T r_{0} )
  double beta = sqrt(dot_product(zz, vv));
  const double beta0 = beta;

  // v_{1} = r_{0} / beta_{0}
  vv *= 1 / beta;
  // pr_{1} = LMP r_{0} / beta_{0}
  pr *= 1 / beta;
  // z_{1} = z_{0} / beta_{0}
  zz *= 1 / beta;

  // hvecs[0] = pr_{1} --> required for solution
  ritzPairs.zVEC().emplace_back(new CtrlInc_(pr));
  // zvecs[0] = z_{1} ---> for re-orthogonalization
  ritzPairs.tVEC().emplace_back(new CtrlInc_(zz));
  // vvecs[0] = v_{1} ---> for re-orthogonalization
  ritzPairs.vVEC().emplace_back(new CtrlInc_(vv));

  double normReduction = 1.0;

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << "DRPLanczos Starting Iteration " << jiter + 1 << std::endl;

    // v_{i+1} = ( pr_{i} + H^T R^{-1} H z_{i} ) - beta * v_{i-1}
    HtRinvH.multiply(zz, vv);
    B.multiplyId(pr, mm);
    vv += mm;
    if (jiter > 0) vv.axpy(-beta, ritzPairs.vVEC(jiter - 1));

    // alpha_{i} = v_{i+1}^T z_{i}
    double alpha = dot_product(zz, vv);

    // v_{i+1} = v_{i+1} - alpha_{i} v_{i}
    vv.axpy(-alpha, ritzPairs.vVEC(jiter));  // vv = vv - alpha * v_j

    // Re-orthogonalization
    for (int jj = 0; jj < jiter; ++jj) {
      double proj = dot_product(vv, ritzPairs.tVEC(jj));
      vv.axpy(-proj, ritzPairs.vVEC(jj));
    }

    // z_{i+1} = B LMP v_{i+1}
    lmp_.multiply(vv, pr);
    B.multiplyB(pr, zz);

    // beta_{i+1} = sqrt( zz_{i+1}^t, vv_{i+1} )
    beta = sqrt(dot_product(zz, vv));

    // v_{i+1} = v_{i+1} / beta_{i+1}
    vv *= 1 / beta;
    // pr_{i+1} = pr_{i+1} / beta_{i+1}
    pr *= 1 / beta;
    // z_{i+1} = z_{i+1} / beta_{i+1}
    zz *= 1 / beta;

    // hvecs[i+1] = pr_{i+1}
    ritzPairs.zVEC().emplace_back(new CtrlInc_(pr));
    // zvecs[i+1] = z_{i+1}
    ritzPairs.tVEC().emplace_back(new CtrlInc_(zz));
    // vvecs[i+1] = v_{i+1}
    ritzPairs.vVEC().emplace_back(new CtrlInc_(vv));

    ritzPairs.alphas().push_back(alpha);

    if (jiter == 0) {
      ss.push_back(beta0 / alpha);
      dd.push_back(beta0);
    } else {
      // Solve the tridiagonal system T_{i} s_{i} = beta0 * e_1
      dd.push_back(beta0 * dot_product(ritzPairs.tVEC(0), vv));
      TriDiagSolve(ritzPairs.alphas(), ritzPairs.betas(), dd, ss);
    }

    ritzPairs.betas().push_back(beta);

    // Compute the quadratic cost function
    // J[du_{i}] = J[0] - 0.5 s_{i}^T Z_{i}^T r_{0}
    // Jb[du_{i}] = 0.5 s_{i}^T H_{i}^T Z_{i} s_{i} + gradJb^T Z_{i} s_{i}
    double costJ = costJ0;

    // Calculate the solution (dxh = Binv dx)
    dx.zero();
    dxh.zero();
    for (unsigned int jj = 0; jj < ss.size(); ++jj) {
      dx.axpy(ss[jj], ritzPairs.tVEC(jj));
      dxh.axpy(ss[jj], ritzPairs.zVEC(jj));
    }
    double costJb =
        costJ0Jb + dot_product(dx, gradJb) + 0.5 * dot_product(dx, dxh);

    costJ -= 0.5 * dot_product(dx, rr);
    double costJoJc = costJ - costJb;

    // Gradient norm in precond metric --> sqrt(r'z) --> beta * s_{i}
    double rznorm = beta * std::abs(ss[jiter]);
    normReduction = rznorm / beta0;

    Log::info() << "DRPLanczos end of iteration " << jiter + 1 << std::endl
                << "  Norm reduction (" << std::setw(2) << jiter + 1
                << ") = " << util::full_precision(normReduction) << std::endl
                << "  Quadratic cost function: J   (" << std::setw(2)
                << jiter + 1 << ") = " << util::full_precision(costJ)
                << std::endl
                << "  Quadratic cost function: Jb  (" << std::setw(2)
                << jiter + 1 << ") = " << util::full_precision(costJb)
                << std::endl
                << "  Quadratic cost function: JoJc(" << std::setw(2)
                << jiter + 1 << ") = " << util::full_precision(costJoJc)
                << std::endl
                << std::endl;

    if (normReduction < tolerance) {
      Log::info() << "DRPLanczos: Achieved required reduction in residual norm."
                  << std::endl;
      break;
    }
  }

  // Calculate the solution (dxh = Binv dx)
  dx.zero();
  dxh.zero();
  for (unsigned int jj = 0; jj < ss.size(); ++jj) {
    dx.axpy(ss[jj], ritzPairs.tVEC(jj));
    dxh.axpy(ss[jj], ritzPairs.zVEC(jj));
  }

  // Process Ritz pairs
  ritzPairs.process(config_, "primal");

  return normReduction;
}

// -----------------------------------------------------------------------------

}  // namespace oops
