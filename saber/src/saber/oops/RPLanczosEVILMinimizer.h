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

#include "eckit/mpi/Comm.h"

#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualMinimizer.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HBHtMatrix.h"
#include "oops/assimilation/RinvMatrix.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/util/Logger.h"
#include "oops/util/formats.h"

#include "saber/oops/RitzPairs.h"

#include "util/dot_product.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class RPLanczosEVILMinimizer : public oops::DualMinimizer<MODEL> {
  using CostFct_ = oops::CostFunction<MODEL>;
  using Dual_ = oops::DualVector<MODEL>;
  using HBHt_ = oops::HBHtMatrix<MODEL>;
  using Rinv_ = oops::RinvMatrix<MODEL>;

 public:
  const std::string classname() const override
    { return "RPLanczosEVILMinimizer"; }
  RPLanczosEVILMinimizer(const eckit::Configuration &conf, const CostFct_ &J)
      : oops::DualMinimizer<MODEL>(J), conf_(conf) {}
  ~RPLanczosEVILMinimizer() {}

 private:
  double solve(Dual_ &, double &, Dual_ &, const HBHt_ &, const Rinv_ &,
               const int &, const double &, Dual_ &, const double &) override;

  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
double RPLanczosEVILMinimizer<MODEL>::solve(Dual_ &vv, double &vvp, Dual_ &rr,
                                        const HBHt_ &HBHt, const Rinv_ &Rinv,
                                        const int &maxiter,
                                        const double &tolerance, Dual_ &dy,
                                        const double &sigma) {
  oops::IdentityMatrix<Dual_> precond;

  Dual_ zz(vv);
  Dual_ ww(vv);
  Dual_ tt(vv);
  Dual_ vold(vv);
  Dual_ v(vv);

  // augmented part of the vectors
  double rrp = 1.0;
  double zzp;
  double wwp;
  double ttp = 0.0;
  double voldp = 0.0;
  double vp = 1.0;

  RitzPairs<Dual_> ritzPairs;
  std::vector<double> vpVEC;
  std::vector<double> tpVEC;
  std::vector<double> zpVEC;
  std::vector<double> y;
  std::vector<double> d;

  // zzaug = Gaug  rraug
  precond.multiply(rr, zz);
  zzp = rrp;

  // ttaug = HBHtaug zzaug
  HBHt.multiply(zz, tt);
  tt.axpy(zzp, dy);
  ttp = dot_product(dy, zz) + sigma * zzp;

  double normReduction = 1.0;
  double beta0 = sqrt(dot_product(rr, tt) + rrp * ttp);
  double beta = 0.0;

  vold.zero();
  v = rr;
  vp = rrp;
  v *= 1 / beta0;
  vp *= 1 / beta0;
  tt *= 1 / beta0;
  zz *= 1 / beta0;
  ttp *= 1 / beta0;
  zzp *= 1 / beta0;

  ritzPairs.vVEC().clear();
  ritzPairs.zVEC().clear();
  ritzPairs.tVEC().clear();
  vpVEC.clear();
  zpVEC.clear();
  tpVEC.clear();

  ritzPairs.betas().clear();
  ritzPairs.alphas().clear();

  ritzPairs.vVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(v)));
  ritzPairs.zVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(zz)));
  ritzPairs.tVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(tt)));
  vpVEC.push_back(vp);
  zpVEC.push_back(zzp);
  tpVEC.push_back(ttp);

  int jiter;
  oops::Log::info() << std::endl;
  for (jiter = 0; jiter < maxiter; ++jiter) {
    oops::Log::info() << "RPLanczos Starting Iteration " << jiter + 1 << std::endl;

    // ww = (RinvHBHt + I) zz - beta * vold
    Rinv.multiply(tt, ww);
    ww += zz;
    wwp = zzp;
    ww.axpy(-beta, vold);
    wwp -= beta * voldp;

    double alpha = dot_product(tt, ww) + ttp * wwp;

    ww.axpy(-alpha, v);  // w = w - alpha * v
    wwp -= alpha * vp;

    // Re-orthogonalization
    for (int iiter = 0; iiter < jiter; ++iiter) {
      double proj = dot_product(ww, ritzPairs.tVEC(iiter)) + wwp * tpVEC[iiter];
      ww.axpy(-proj, ritzPairs.vVEC(iiter));
      wwp -= proj * vpVEC[iiter];
    }

    // wwaug = Gaug  zzaug
    precond.multiply(ww, zz);
    zzp = wwp;

    // ttaug = HBHtaug zzaug
    HBHt.multiply(zz, tt);
    tt.axpy(zzp, dy);
    ttp = dot_product(dy, zz) + sigma * zzp;

    beta = sqrt(dot_product(tt, ww) + ttp * wwp);

    vold = v;
    voldp = vp;
    v = ww;
    vp = wwp;
    v *= 1 / beta;
    vp *= 1 / beta;
    tt *= 1 / beta;
    zz *= 1 / beta;
    ttp *= 1 / beta;
    zzp *= 1 / beta;

    ritzPairs.vVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(v)));
    ritzPairs.zVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(zz)));
    ritzPairs.tVEC().emplace_back(std::unique_ptr<Dual_>(new Dual_(tt)));
    vpVEC.push_back(vp);
    zpVEC.push_back(zzp);
    tpVEC.push_back(ttp);

    ritzPairs.alphas().push_back(alpha);

    y.clear();
    if (jiter == 0) {
      y.push_back(beta0 / alpha);
    } else {
      // Solve the tridiagonal system T_jiter y_jiter = beta0 * e_1
      d.clear();
      for (int iiter = 0; iiter <= jiter; ++iiter) {
        d.push_back(beta0 *
                    (dot_product(ritzPairs.tVEC(0), ritzPairs.vVEC(iiter)) +
                     tpVEC[0] * vpVEC[iiter]));
      }
      oops::TriDiagSolve(ritzPairs.alphas(), ritzPairs.betas(), d, y);
    }

    // Gradient norm in precond metric --> sqrt(r't) --> beta * y(jiter)
    double rznorm = beta * std::abs(y[jiter]);

    normReduction = rznorm / beta0;

    ritzPairs.betas().push_back(beta);

    oops::Log::info() << "RPLanczos end of iteration " << jiter + 1
                << ". Norm reduction= " << util::full_precision(normReduction)
                << std::endl
                << std::endl;

    if (normReduction < tolerance) {
      oops::Log::info() << "RPLanczos: Achieved required reduction in residual norm."
                  << std::endl;
      break;
    }
  }

  // Calculate the solution (xh = Binv x)
  for (int iiter = 0; iiter < jiter; ++iiter) {
    vv.axpy(y[iiter], ritzPairs.zVEC(iiter));
    vvp += y[iiter] * zpVEC[iiter];
  }

  // Process Ritz pairs
  ritzPairs.process(conf_);

  oops::Log::info() << "RPLanczos: end" << std::endl;

  return normReduction;
}

// -----------------------------------------------------------------------------

}  // namespace saber
