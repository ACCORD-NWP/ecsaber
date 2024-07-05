/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/TriDiagSpectrum.h"

#include "oops/util/ConfigFunctions.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename VECTOR>
class RitzPairs {
 public:
  // Constructor and destructor
  RitzPairs(){};
  ~RitzPairs() {}

  // Data accessors
  const std::vector<std::shared_ptr<VECTOR>> &vVEC() const { return vVEC_; }
  std::vector<std::shared_ptr<VECTOR>> &vVEC() { return vVEC_; }
  const VECTOR &vVEC(const int &ii) const { return *vVEC_[ii]; }
  VECTOR &vVEC(const int &ii) { return *vVEC_[ii]; }
  const std::vector<std::shared_ptr<VECTOR>> &zVEC() const { return zVEC_; }
  std::vector<std::shared_ptr<VECTOR>> &zVEC() { return zVEC_; }
  const VECTOR &zVEC(const int &ii) const { return *zVEC_[ii]; }
  VECTOR &zVEC(const int &ii) { return *zVEC_[ii]; }
  const std::vector<std::shared_ptr<VECTOR>> &tVEC() const { return tVEC_; }
  std::vector<std::shared_ptr<VECTOR>> &tVEC() { return tVEC_; }
  const VECTOR &tVEC(const int &ii) const { return *tVEC_[ii]; }
  VECTOR &tVEC(const int &ii) { return *tVEC_[ii]; }

  const std::vector<double> &alphas() const { return alphas_; }
  std::vector<double> &alphas() { return alphas_; }
  const std::vector<double> &betas() const { return betas_; }
  std::vector<double> &betas() { return betas_; }

  // Process Ritz pairs
  void process(const eckit::Configuration &);

 private:
  std::vector<std::shared_ptr<VECTOR>> vVEC_;
  std::vector<std::shared_ptr<VECTOR>> zVEC_;
  std::vector<std::shared_ptr<VECTOR>> tVEC_;
  std::vector<double> alphas_;
  std::vector<double> betas_;
};

// -----------------------------------------------------------------------------

template <typename VECTOR>
void RitzPairs<VECTOR>::process(const eckit::Configuration &conf) {
  oops::Log::trace() << "RitzPairs::process start" << std::endl;

  // Get number of Ritz pairs
  const int jiter = vVEC_.size() - 1;

  // Get Ritz convergence tolerance
  double ritzTol = conf.getDouble("ritz tolerance", 1.0);

  // Compute eigenpairs
  oops::Log::info() << "Compute eigenpairs" << std::endl;
  std::vector<double> evals;
  std::vector<std::vector<double>> evecs;
  oops::TriDiagSpectrum(alphas_, betas_, evals, evecs);

  // Compute valid Ritz pairs
  oops::Log::info() << "Compute valid Ritz pairs" << std::endl;
  eckit::LocalConfiguration ritzConfTemplate(conf, "ritz vectors");
  std::vector<double> ritzVal;
  for (int iiter1 = jiter - 1; iiter1 >= 0; --iiter1) {
    // Compute convergence bound
    double ebnd =
        std::abs(betas_[iiter1] * evecs[jiter - 1][iiter1]) / alphas_[iiter1];

    if (ebnd < ritzTol) {
      oops::Log::test() << "Accept Ritz pair #" << (jiter - iiter1) << std::endl;
      oops::Log::info() << "Eigen value: " << evals[iiter1] << " (" << ebnd << " < "
                        << ritzTol << ")" << std::endl;

      // Save Ritz value
      ritzVal.push_back(evals[iiter1]);

      // Compute first Ritz vector
      VECTOR ritzBar(*vVEC_[0], false);
      ritzBar.zero();
      for (int iiter2 = 0; iiter2 < jiter; ++iiter2) {
        VECTOR tmpBar(*vVEC_[iiter2]);
        tmpBar *= evecs[iiter1][iiter2];
        ritzBar += tmpBar;
      }

      // Write first Ritz vector
      eckit::LocalConfiguration ritzBarConf(ritzConfTemplate);
      util::seekAndReplace(ritzBarConf, "%pattern%", "bar");
      util::seekAndReplace(ritzBarConf, "%iteration%", jiter - iiter1, 0);
      ritzBar.write(ritzBarConf);

      // Optionally deal with a second vector
      if (tVEC_.size() > 0 || zVEC_.size() > 0) {
        // Compute second Ritz vector
        VECTOR ritzHat(*vVEC_[0], false);
        ritzHat.zero();
        for (int iiter2 = 0; iiter2 < jiter; ++iiter2) {
          // Second vector
          VECTOR tmpHat(*vVEC_[iiter2], false);
          if (tVEC_.size() > 0) {
            tmpHat = *tVEC_[iiter2];
          } else {
            tmpHat = *zVEC_[iiter2];
          }
          tmpHat *= evecs[iiter1][iiter2];
          ritzHat += tmpHat;
        }

        // Write second Ritz vector
        eckit::LocalConfiguration ritzHatConf(ritzConfTemplate);
        util::seekAndReplace(ritzHatConf, "%pattern%", "hat");
        util::seekAndReplace(ritzHatConf, "%iteration%", jiter - iiter1, 0);
        ritzHat.write(ritzHatConf);
      }
    } else {
      oops::Log::test() << "Dismiss Ritz pair #" << (jiter - iiter1) << std::endl;
      oops::Log::info() << "Eigen value: " << evals[iiter1] << " (" << ebnd << " < "
                        << ritzTol << ")" << std::endl;
    }
  }

  // Write Ritz values
  const eckit::mpi::Comm &comm(eckit::mpi::comm());
  if (comm.rank() == 0) {
    std::ofstream ritzValFile(conf.getString("eigenvalues").c_str());
    if (!ritzValFile.is_open()) ABORT("Unable to open eigenvalues file");
    for (const auto &val : ritzVal) {
      ritzValFile << val << std::endl;
    }
    ritzValFile.close();
  }

  oops::Log::trace() << "RitzPairs::process done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
