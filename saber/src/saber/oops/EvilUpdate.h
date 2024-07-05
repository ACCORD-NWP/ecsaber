/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <fstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/ControlVector.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/assimilation/RinvMatrix.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/Departures.h"
#include "oops/base/Ensemble.h"
#include "oops/base/Increment4D.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/Variables.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/abor1_cpp.h"

namespace saber {

template <typename MODEL>
class EvilUpdate : public oops::Application {
  using BMatrix_ = oops::BMatrix<MODEL>;
  using CtrlInc_ = oops::ControlIncrement<MODEL>;
  using CtrlVariable_ = oops::ControlVariable<MODEL>;
  using CtrlVec_ = oops::ControlVector<MODEL>;
  using Departures_ = oops::Departures<MODEL>;
  using Dual_ = oops::DualVector<MODEL>;
  using Ensemble_ = oops::Ensemble<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using HMatrix_ = oops::HMatrix<MODEL>;
  using HtMatrix_ = oops::HtMatrix<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using RinvMatrix_ = oops::RinvMatrix<MODEL>;
  using State_ = oops::State<MODEL>;
  using Variables_ = oops::Variables<MODEL>;

 public:
  // -----------------------------------------------------------------------------
  EvilUpdate() {
    oops::instantiateCostFactory<MODEL>();
    oops::instantiateCovarFactory<MODEL>();
    oops::instantiateMinFactory<MODEL>();
    oops::instantiateObsErrorFactory<MODEL>();
    oops::instantiateTlmFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~EvilUpdate() {}
  // -----------------------------------------------------------------------------
  /*
Implementation of EVIL algorithms to update ensemble members (similar function
to an EDA).

A minimization is performed in control (IncrCtlVec), primal (Increment) or dual
(Observations) space, using either the SQRTPLanczos, PLanczos or the RPLanczos
minimizer, respectively. Ritz pairs are computed as a by-product, each Ritz pair
containing:
- an eigenvalue \theta
- one or two eigenvectors:
  + z^c in control space,
  + (\hat{z}^p, \bar{z}^p) in primal space,
  + (\hat{z}^d, \bar{z}^d) in dual space

The eigen vectors are linked:
- In primal space: \hat{z}^p = B \bar{z}^p
- In dual space: \hat{z}^d = HBH^T \bar{z}^d
- Between control and primal spaces: \hat{z}^p = U z^c and z^c = U^T \bar{z}^p
- Between primal and dual spaces: \hat{z}^d = H \hat{z}^p and \hat{z}^p = BH^T
\bar{z}^d and \bar{z}^p = H^T \bar{z}^d

EVIL updates for 1 <= k <= N:
 - S-EVIL: x^a_k = x^b_k + \sum_{l=1}^{N_e} w^S_{kl} \hat{z}^p_l
 - D-EVIL: x^a_k = x^b_k + \sum_{l=1}^{N_e} w^D_{kl} \hat{z}^p_l
 - R-EVIL: x^a_k = U \ksi_k + \sum_{l=1}^{N_e} w^R_{kl} \hat{z}^p_l

where:
- N is the number of ensemble members (randomized or read)
- N_e is the number of pre-computed Ritz pairs
- \hat{z}^p_l is the lth first Ritz vector in primal space
- \ksi_k is a random Gaussian vector of zero mean and unit variance
- The weights are given by:
  + S-EVIL: w^S_{kl} = 1/\theta_l (\hat{z}^d_l)^T R{-1} (y_k - Hx^b_k)
    where y_k = R^{1/2} \eta is a random observations perturbation
  + D-EVIL: w^D_{kl} = -(1 - 1/\sqrt{theta_l}) (\bar{z}^p_l)^T x^b_k
    which is not available in control space since \bar{z}^p_l would require
B^{-1}
  + S-EVIL: w^S_{kl} = -(1 - 1/\sqrt{theta_l}) (z^c_l)^T \ksi_k

Reference:
https://journals.ametsoc.org/view/journals/mwre/144/10/mwr-d-15-0252.1.xml
  */
  int execute(const eckit::Configuration &fullConfig) const {
    // Get EVIL parameters
    const eckit::LocalConfiguration evilConfig(fullConfig, "evil");

    // Setup EVIL update type
    const std::string filter = evilConfig.getString("filter");
    ASSERT(filter == "S" || filter == "D" || filter == "R");

    // Setup EVIL space ("control, "primal" or "dual")
    const std::string solverSpace = evilConfig.getString("solver space");
    ASSERT(solverSpace == "control" || solverSpace == "primal" ||
           solverSpace == "dual");

    // Check consistency (D-EVIL not available in control space)
    ASSERT(!(filter == "D" && solverSpace == "control"));

    // Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    // Setup model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    // Setup variables
    const eckit::LocalConfiguration varConfig(fullConfig, "variables");
    const Variables_ vars(varConfig);

    // Setup cost function
    const eckit::LocalConfiguration cfConf(fullConfig, "cost_function");
    std::unique_ptr<oops::CostFunction<MODEL> > J(
        oops::CostFactory<MODEL>::create(cfConf, resol, model));
    oops::Log::trace() << "EvilUpdate: cost function has been set up" << std::endl;

    // Initialize first guess from background
    const CtrlVariable_ xx(J->jb().getBackground());

    // Linearize cost function
    oops::PostProcessor<State_> post;
    std::vector<eckit::LocalConfiguration> iterconfs;
    fullConfig.get("variational.iteration", iterconfs);
    J->linearize(xx, iterconfs[0], post);

    // Setup variance reduction
    const eckit::LocalConfiguration varRedConfTemplate(fullConfig,
                                                       "variance reduction");
    std::vector<double> varianceNorm;

    // Define linearizing state
    const State_ xb(xx.state()[0]);

    // Setup H matrix
    std::unique_ptr<HMatrix_> HMatrix;
    if (filter == "S" || solverSpace == "dual") {
      oops::Log::info() << "Setup H matrix" << std::endl;
      HMatrix.reset(new HMatrix_(*J));
    }

    // Setup Ht matrix
    std::unique_ptr<HtMatrix_> HtMatrix;
    if (solverSpace == "dual") {
      oops::Log::info() << "Setup Ht matrix" << std::endl;
      HtMatrix.reset(new HtMatrix_(*J));
    }

    // Setup R matrix
    std::unique_ptr<RinvMatrix_> RinvMatrix;
    if (filter == "S") {
      oops::Log::info() << "Setup R matrix" << std::endl;
      RinvMatrix.reset(new RinvMatrix_(*J));
    }

    // Setup B matrix
    std::unique_ptr<BMatrix_> BMatrix;
    if (filter == "R" || solverSpace == "control" || solverSpace == "dual") {
      oops::Log::info() << "Setup B matrix" << std::endl;
      BMatrix.reset(new BMatrix_(*J));
    }

    // Read analysis if available
    std::unique_ptr<State_> xa;
    if (fullConfig.has("analysis")) {
      oops::Log::info() << "Read analysis" << std::endl;
      eckit::LocalConfiguration anaConf(fullConfig, "analysis");
      xa.reset(new State_(resol, model, anaConf));
    }

    // Setup ensemble of backgrounds perturbations
    const eckit::LocalConfiguration ensConfig(fullConfig,
                                              "ensemble of backgrounds");
    Ensemble_ Xb(xb.validTime(), ensConfig);
    size_t nens = Xb.size();
    std::vector<CtrlVec_> ksi;
    if (filter == "R") {
      // Generate random control vectors
      oops::Log::info() << "Generate random control vectors" << std::endl;
      for (size_t ie = 0; ie < nens; ++ie) {
        ksi.push_back(CtrlVec_(J->jb()));
        ksi[ie].zero();
        ksi[ie].state().random();
      }

      // Compute random control vectors mean
      oops::Log::info() << "Compute random control vectors mean" << std::endl;
      CtrlVec_ ksiMean(J->jb());
      ksiMean.zero();
      for (size_t ie = 0; ie < nens; ++ie) {
        ksiMean += ksi[ie];
      }
      const double meanNorm = 1.0 / static_cast<double>(nens);
      ksiMean *= meanNorm;

      // Remove random control vectors mean
      oops::Log::info() << "Remove random control vectors mean" << std::endl;
      for (size_t ie = 0; ie < nens; ++ie) {
        ksi[ie] -= ksiMean;
      }

      // Apply background error covariance square-root to control vectors
      oops::Log::info()
          << "Apply background error covariance square-root to control vectors"
          << std::endl;
      Xb.build(xb, resol);
      CtrlInc_ dx(J->jb());
      for (size_t ie = 0; ie < nens; ++ie) {
        BMatrix->multiplySqrt(ksi[ie], dx);
        Xb[ie] = dx.state()[dx.state().first()];
      }
    } else {
      // Read ensemble of backgrounds perturbations
      oops::Log::info() << "Read ensemble of backgrounds" << std::endl;
      Xb.linearize(xb, resol);
    }

    // Compute background variance
    Increment_ variance = Xb.variance();
    eckit::LocalConfiguration varRedConf(varRedConfTemplate);
    util::seekAndReplace(varRedConf, "%iteration%", 0);
    variance.write(varRedConf);
    varianceNorm.push_back(variance.norm());

    // Convert ensemble of backgrounds perturbations from primal to dual space
    std::vector<Dual_> HXb;
    if (filter == "S" || solverSpace == "dual") {
      // Setup ControlIncrement and DualVector
      CtrlInc_ dx(J->jb());
      Dual_ dy;
      for (unsigned jj = 0; jj < J->nterms(); ++jj) {
        dy.append(J->jterm(jj).newDualVector());
      }

      for (size_t ie = 0; ie < nens; ++ie) {
        // Apply H matrix to Xb members
        dx.state()[dx.state().first()] = Xb[ie];
        HMatrix->multiply(dx, dy);
        HXb.push_back(dy);
      }
    }

    // Compute difference between observations and background perturbations
    std::vector<Dual_> ypertVec;
    if (filter == "S") {
      // Generate and store normalized observation perturbations
      for (size_t ie = 0; ie < nens; ++ie) {
        Dual_ ypert;
        for (unsigned jj = 0; jj < J->nterms(); ++jj) {
          if (J->jterm(jj).classname() == "Jo") {
            ypert.append(J->jterm(jj).randomizeCovar());
          } else {
            ypert.append(J->jterm(jj).newDualVector());
          }
        }
        ypertVec.push_back(ypert);
      }

      // Compute and subtract observation perturbations mean
      Dual_ ypertMean;
      for (unsigned jj = 0; jj < J->nterms(); ++jj) {
        ypertMean.append(J->jterm(jj).newDualVector());
      }
      ypertMean.zero();
      for (size_t ie = 0; ie < nens; ++ie) {
        ypertMean += ypertVec[ie];
      }
      ypertMean *= 1.0 / static_cast<double>(nens);
      for (size_t ie = 0; ie < nens; ++ie) {
        ypertVec[ie] -= ypertMean;
      }

      for (size_t ie = 0; ie < nens; ++ie) {
        // Normalize perturbations
        ypertVec[ie] *= 1.0 / std::sqrt(static_cast<double>(nens - 1));

        // Compute perturbations difference
        ypertVec[ie] -= HXb[ie];

        // Apply R inverse on difference
        RinvMatrix->multiply(ypertVec[ie], HXb[ie]);
      }
    }

    // Read eigenvalues
    oops::Log::info() << "Read eigenvalues" << std::endl;
    std::vector<double> evals;
    size_t niterMax;
    const eckit::mpi::Comm &comm(eckit::mpi::comm());
    if (comm.rank() == 0) {
      std::string line;
      std::ifstream evalsFile(evilConfig.getString("eigenvalues").c_str());
      if (!evalsFile.is_open()) ABORT("Unable to open eigenvalues file");
      while (evalsFile.peek() != EOF) {
        std::getline(evalsFile, line);
        evals.push_back(std::stod(line));
      }
      evalsFile.close();
      niterMax = evals.size();
      comm.broadcast(niterMax, 0);
      comm.broadcast(evals, 0);
    } else {
      comm.broadcast(niterMax, 0);
      evals.resize(niterMax);
      comm.broadcast(evals, 0);
    }

    // Get number of Ritz pairs
    const size_t niter = evilConfig.getInt("number of ritz pairs", niterMax);
    ASSERT(niter <= niterMax);
    oops::Log::info() << "Number of Ritz pairs: " << niter << " (max. " << niterMax
                      << ")" << std::endl;

    // Copy ensemble of backgrounds to initialize ensemble of analyses
    Ensemble_ Xa(Xb);

    // Loop over Ritz pairs
    eckit::LocalConfiguration ritzConfTemplate(evilConfig, "ritz vectors");
    for (size_t iiter = 0; iiter < niter; ++iiter) {
      // Prepare Ritz vectors configurations
      eckit::LocalConfiguration ritzCtlConf(ritzConfTemplate);
      util::seekAndReplace(ritzCtlConf, "%iteration%", iiter+1, 0);
      eckit::LocalConfiguration ritzHatConf(ritzConfTemplate);
      util::seekAndReplace(ritzHatConf, "%pattern%", "hat");
      util::seekAndReplace(ritzHatConf, "%iteration%", iiter+1, 0);
      eckit::LocalConfiguration ritzBarConf(ritzConfTemplate);
      util::seekAndReplace(ritzBarConf, "%pattern%", "bar");
      util::seekAndReplace(ritzBarConf, "%iteration%", iiter+1, 0);

      // Initialize Ritz vectors
      std::unique_ptr<CtrlVec_> ritzCtl;
      std::unique_ptr<CtrlInc_> ritzHatIncr;
      std::unique_ptr<CtrlInc_> ritzBarIncr;
      std::unique_ptr<Dual_> ritzHatDual;
      std::unique_ptr<Dual_> ritzBarDual;

      if (solverSpace == "control") {
        // Read Ritz vector in control space
        ritzCtl.reset(new CtrlVec_(J->jb()));
        ritzCtl->read(ritzCtlConf);

        // Compute first Ritz vector in primal space
        ritzHatIncr.reset(new CtrlInc_(J->jb()));
        BMatrix->multiplySqrt(*ritzCtl, *ritzHatIncr);

        if (filter == "S") {
          // Compute first Ritz vector in dual space
          ritzHatDual.reset(new Dual_());
          for (unsigned jj = 0; jj < J->nterms(); ++jj) {
            ritzHatDual->append(J->jterm(jj).newDualVector());
          }
          HMatrix->multiply(*ritzHatIncr, *ritzHatDual);
        }
      } else if (solverSpace == "primal") {
        // Read first Ritz vector in primal space
        ritzHatIncr.reset(new CtrlInc_(J->jb()));
        Increment_ incr(resol, vars, xb.validTime());
        incr.read(ritzHatConf);
        ritzHatIncr->state()[ritzHatIncr->state().first()] = incr;

        if (filter == "D" || filter == "R") {
          // Read second Ritz vector in primal space
          ritzBarIncr.reset(new CtrlInc_(J->jb()));
          incr.read(ritzBarConf);
          ritzBarIncr->state()[ritzBarIncr->state().first()] = incr;
        }

        if (filter == "S") {
          // Compute first vector in dual space
          ritzHatDual.reset(new Dual_());
          for (unsigned jj = 0; jj < J->nterms(); ++jj) {
            ritzHatDual->append(J->jterm(jj).newDualVector());
          }
          HMatrix->multiply(*ritzHatIncr, *ritzHatDual);
        }

        if (filter == "R") {
          // Compute Ritz vector in control space
          ritzCtl.reset(new CtrlVec_(J->jb()));
          BMatrix->multiplySqrtTrans(*ritzBarIncr, *ritzCtl);
        }
      } else if (solverSpace == "dual") {
        if (filter == "S") {
          // Read first Ritz vector in primal space
          ritzHatDual.reset(new Dual_());
          for (unsigned jj = 0; jj < J->nterms(); ++jj) {
            ritzHatDual->append(J->jterm(jj).newDualVector());
          }
          ritzHatDual->read(ritzHatConf);
        }

        // Read second Ritz vector in primal space
        ritzBarDual.reset(new Dual_());
        for (unsigned jj = 0; jj < J->nterms(); ++jj) {
          ritzBarDual->append(J->jterm(jj).newDualVector());
        }
        ritzBarDual->read(ritzBarConf);

        // Compute Ritz vectors in primal space
        ritzBarIncr.reset(new CtrlInc_(J->jb()));
        HtMatrix->multiply(*ritzBarDual, *ritzBarIncr);
        ritzHatIncr.reset(new CtrlInc_(J->jb()));
        BMatrix->multiply(*ritzBarIncr, *ritzHatIncr);

        if (filter == "R") {
          // Compute Ritz vector in control space
          ritzCtl.reset(new CtrlVec_(J->jb()));
          BMatrix->multiplySqrtTrans(*ritzBarIncr, *ritzCtl);
        }
      }

      // Update analysis perturbations
      for (size_t ie = 0; ie < nens; ++ie) {
        // Compute weight
        double weight;
        if (filter == "S") {
          weight =
              (1.0 / evals[iiter]) * ritzHatDual->dot_product_with(HXb[ie]);
        } else if (filter == "D") {
          weight = -(1.0 - 1.0 / std::sqrt(evals[iiter])) *
                   ritzBarIncr->state()[ritzBarIncr->state().first()]
                       .dot_product_with(Xb[ie]);
        } else if (filter == "R") {
          weight = -(1.0 - 1.0 / std::sqrt(evals[iiter])) *
                   ritzCtl->dot_product_with(ksi[ie]);
        }

        // Update analysis perturbation
        Xa[ie].axpy(weight, ritzHatIncr->state()[ritzHatIncr->state().first()]);
      }

      // Compute analysis variance
      Increment_ variance = Xa.variance();
      eckit::LocalConfiguration varRedConf(varRedConfTemplate);
      util::seekAndReplace(varRedConf, "%iteration%", iiter, 0);
      variance.write(varRedConf);
      varianceNorm.push_back(variance.norm());

      // Compute variance reduction
      double varianceReduction =
          varianceNorm[iiter + 1] / varianceNorm[0] * 100.0;
      std::streamsize ss = oops::Log::test().precision();
      oops::Log::test() << "Update with Ritz pair #" << iiter << ": variance is "
                        << std::fixed << std::setprecision(2) << varianceReduction
                        << "% of the initial variance" << std::endl;
      oops::Log::test().precision(ss);
    }

    // Compute and write ensemble of analyses
    oops::Log::info() << "Compute and write ensemble of analyses" << std::endl;
    eckit::LocalConfiguration xaPertConfTemplate(fullConfig,
                                                 "ensemble of analyses");
    for (size_t ie = 0; ie < nens; ++ie) {
      Xa[ie] *= std::sqrt(static_cast<double>(nens - 1));
      eckit::LocalConfiguration xaPertConf(xaPertConfTemplate);
      util::seekAndReplace(xaPertConf, "%member%", ie+1, 0);
      if (fullConfig.has("analysis")) {
        // Add analysis perturbations to the analysis and write
        State_ xaPert(*xa);
        xaPert += Xa[ie];
        xaPert.write(xaPertConf);
      } else {
        // Write analysis perturbations directly
        Xa[ie].write(xaPertConf);
      }
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EvilUpdate<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace saber
