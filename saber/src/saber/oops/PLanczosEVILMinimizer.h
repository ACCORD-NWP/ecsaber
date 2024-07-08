/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/HessianMatrix.h"
#include "oops/assimilation/PrimalMinimizer.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/util/Logger.h"
#include "oops/util/formats.h"

#include "saber/oops/RitzPairs.h"

#include "util/dot_product.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class PLanczosEVILMinimizer : public oops::PrimalMinimizer<MODEL> {
  using Bmat_ = oops::BMatrix<MODEL>;
  using CostFct_ = oops::CostFunction<MODEL>;
  using CtrlInc_ = oops::ControlIncrement<MODEL>;
  using Hessian_ = oops::HessianMatrix<MODEL>;

 public:
  const std::string classname() const override
    { return "PLanczosEVILMinimizer"; }
  PLanczosEVILMinimizer(const eckit::Configuration &conf,
                        const CostFct_ &J)
    : oops::PrimalMinimizer<MODEL>(J), conf_(conf) {}
  ~PLanczosEVILMinimizer() {}

 private:
  double solve(CtrlInc_ &,
               const CtrlInc_ &,
               const Hessian_ &,
               const Bmat_ &,
               const int,
               const double) override;

  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
double PLanczosEVILMinimizer<MODEL>::solve(CtrlInc_ &dx,
                                           const CtrlInc_ &rhs,
                                           const Hessian_ &hessian,
                                           const Bmat_ &B,
                                           const int ninner,
                                           const double gnreduc) {
  // Solve the linear system
  CtrlInc_ zz(dx);
  CtrlInc_ ww(dx);

  RitzPairs<CtrlInc_> ritzPairs;
  std::vector<double> dd;
  std::vector<double> yy;

  // Initial residual r = b - Ax
  CtrlInc_ rr(rhs);
  double xnrm2 = dot_product(dx, dx);
  if (xnrm2 != 0) {
    hessian.multiply(dx, zz);
    rr -= zz;
  }

  // z = B r
  B.multiply(rr, zz);

  double reduc = 1.0;
  double beta0 = sqrt(dot_product(rr, zz));
  double beta = 0.0;

  CtrlInc_ vv(rr);
  vv *= 1 / beta0;
  zz *= 1 / beta0;

  // zVEC[0] = z_1 ---> required for re-orthogonalization
  ritzPairs.zVEC().emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));
  // vVEC[0] = v_1 ---> required for re-orthogonalization
  ritzPairs.vVEC().emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));

  int jiter;
  oops::Log::info() << std::endl;
  for (jiter = 0; jiter < ninner; ++jiter) {
    oops::Log::info() << " PLanczosEVILMinimizer Starting Iteration " << jiter + 1 << std::endl;

    // w = hessian z - beta * vold
    hessian.multiply(zz, ww);  // w = hessian z
    if (jiter > 0) ww.axpy(-beta, ritzPairs.vVEC(jiter - 1));

    double alpha = dot_product(zz, ww);

    ww.axpy(-alpha, vv);  // w = w - alpha * v

    // Re-orthogonalization
    for (int iiter = 0; iiter < jiter; ++iiter) {
      double proj = dot_product(ww, ritzPairs.zVEC(iiter));
      ww.axpy(-proj, ritzPairs.vVEC(iiter));
    }

    B.multiply(ww, zz);  // z = B w

    beta = sqrt(dot_product(zz, ww));

    vv = ww;
    vv *= 1 / beta;
    zz *= 1 / beta;

    // zVEC[jiter+1] = z_jiter
    ritzPairs.zVEC().emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));
    // vVEC[jiter+1] = v_jiter
    ritzPairs.vVEC().emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));

    ritzPairs.alphas().push_back(alpha);

    if (jiter == 0) {
      yy.push_back(beta0 / alpha);
      dd.push_back(beta0);
    } else {
      // Solve the tridiagonal system T_jiter y_jiter = beta0 * e_1
      dd.push_back(beta0 * dot_product(ritzPairs.zVEC(0), vv));
      oops::TriDiagSolve(ritzPairs.alphas(), ritzPairs.betas(), dd, yy);
    }

    // Gradient norm in B metric --> sqrt(r'z) --> beta * y(jiter)
    double rznorm = beta * std::abs(yy[jiter]);

    reduc = rznorm / beta0;

    ritzPairs.betas().push_back(beta);

    oops::Log::info() << "PLanczosEVILMinimizer end of iteration " << jiter + 1
                      << ". Norm reduction= " << util::full_precision(reduc)
                      << std::endl
                      << std::endl;

    if (reduc < gnreduc) {
      oops::Log::info() << "PLanczosEVILMinimizer: Achieved required reduction in residual norm."
                        << std::endl;
      break;
    }
  }

  // Calculate the solution (xh = Binv x)
  for (int iiter = 0; iiter < jiter; ++iiter) {
    dx.axpy(yy[iiter], ritzPairs.zVEC(iiter));
  }

  // Process Ritz pairs
  ritzPairs.process(conf_);

  return reduc;
}

// -----------------------------------------------------------------------------

}  // namespace saber
