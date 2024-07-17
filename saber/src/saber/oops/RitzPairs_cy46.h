/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/TriDiagSpectrum.h"

#include "oops/util/ConfigFunctions.h"
#include "oops/util/formats.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename VECTOR>
class RitzPairs {
 public:
  // Constructor and destructor
  RitzPairs(){};
  ~RitzPairs() {}

  // Data accessors
  const boost::ptr_vector<VECTOR> &vVEC() const { return vVEC_; }
  boost::ptr_vector<VECTOR> &vVEC() { return vVEC_; }
  const VECTOR &vVEC(const int &ii) const { return vVEC_[ii]; }
  VECTOR &vVEC(const int &ii) { return vVEC_[ii]; }
  const boost::ptr_vector<VECTOR> &zVEC() const { return zVEC_; }
  boost::ptr_vector<VECTOR> &zVEC() { return zVEC_; }
  const VECTOR &zVEC(const int &ii) const { return zVEC_[ii]; }
  VECTOR &zVEC(const int &ii) { return zVEC_[ii]; }
  const boost::ptr_vector<VECTOR> &tVEC() const { return tVEC_; }
  boost::ptr_vector<VECTOR> &tVEC() { return tVEC_; }
  const VECTOR &tVEC(const int &ii) const { return tVEC_[ii]; }
  VECTOR &tVEC(const int &ii) { return tVEC_[ii]; }

  const std::vector<double> &alphas() const { return alphas_; }
  std::vector<double> &alphas() { return alphas_; }
  const std::vector<double> &betas() const { return betas_; }
  std::vector<double> &betas() { return betas_; }

  // Process Ritz pairs
  void process(const eckit::Configuration &,
               const std::string &);

 private:
  boost::ptr_vector<VECTOR> vVEC_;
  boost::ptr_vector<VECTOR> zVEC_;
  boost::ptr_vector<VECTOR> tVEC_;
  std::vector<double> alphas_;
  std::vector<double> betas_;
};

// -----------------------------------------------------------------------------

template <typename VECTOR>
void RitzPairs<VECTOR>::process(const eckit::Configuration & conf,
                                const std::string & solverSpace) {
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

      // Compute control or first Ritz vector
      VECTOR ritzBar(vVEC_[0], false);
      ritzBar.zero();
      for (int iiter2 = 0; iiter2 < jiter; ++iiter2) {
        VECTOR tmpBar(vVEC_[iiter2]);
        tmpBar *= evecs[iiter1][iiter2];
        ritzBar += tmpBar;
      }

      // Write control or first Ritz vector
      eckit::LocalConfiguration ritzBarConf(ritzConfTemplate);
      util::seekAndReplace(ritzBarConf, "%pattern%", "bar");
      util::seekAndReplace(ritzBarConf, "%iteration%", jiter - iiter1, 0);
      setMPI(ritzBarConf, eckit::mpi::comm().size());
      ritzBar.state()[ritzBar.state().first()].write(ritzBarConf);

      // Compute second Ritz vector
      VECTOR ritzHat(vVEC_[0], false);
      ritzHat.zero();
      for (int iiter2 = 0; iiter2 < jiter; ++iiter2) {
        // Second vector
        VECTOR tmpHat(vVEC_[iiter2], false);
        tmpHat = tVEC_[iiter2];
        tmpHat *= evecs[iiter1][iiter2];
        ritzHat += tmpHat;
      }

      // Write second Ritz vector
      eckit::LocalConfiguration ritzHatConf(ritzConfTemplate);
      util::seekAndReplace(ritzHatConf, "%pattern%", "hat");
      util::seekAndReplace(ritzHatConf, "%iteration%", jiter - iiter1, 0);
      setMPI(ritzHatConf, eckit::mpi::comm().size());
      ritzHat.state()[ritzHat.state().first()].write(ritzHatConf);
    } else {
      oops::Log::test() << "Dismiss Ritz pair #" << (jiter - iiter1) << std::endl;
      oops::Log::info() << "Eigen value: " << evals[iiter1] << " (" << ebnd << " < "
                        << ritzTol << ")" << std::endl;
    }
  }

  // Write Ritz eigenvalues
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
