/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ObsError.h"

#include "eckit/config/Configuration.h"

#include "util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

static oops::ObsErrorMaker<Traits, ObsError> makerObsErrorDefault_("default");

// -----------------------------------------------------------------------------

ObsError::ObsError(const ObsSpace & obsSpace,
                   const eckit::Configuration & config)
  : lvarqc_(false), cleft_(0.0), cright_(0.0), stddev_(), inverseVariance_(), wghtsqrt_() {
  oops::Log::trace() << classname() << "::ObsError starting" << std::endl;

/// Setup R
  stddev_.reset(new ObsVector(obsSpace));
  const std::string col = config.getString("obserror");
  stddev_->read(col);

  inverseVariance_.reset(new ObsVector(*stddev_));
  *inverseVariance_ *= *inverseVariance_;
  inverseVariance_->invert();

/// Setup W^1/2
  if (config.has("varqc")) {
    lvarqc_ = true;
    cleft_ = config.getDouble("varqc.cleft");
    cright_ = config.getDouble("varqc.cright");
    wghtsqrt_.reset(new ObsVector(obsSpace));
  }

  oops::Log::trace() << classname() << "::ObsError done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsError::setupWeights(const ObsVector & dy) {
  oops::Log::trace() << classname() << "::setupWeights starting" << std::endl;

  if (lvarqc_) {
    // Setup W^1/2
    for (size_t jo = 0; jo < dy.size(); ++jo) {
      // Initialization
      (*wghtsqrt_)(jo) = 1.0;

      if (dy(jo) < cleft_ || dy(jo) > cright_) {
        // rho(x) = x^2/sigma^2         if |x|<=c
        // rho(x) = (2c|x|-c^2)/sigma^2 if |x|>c
        // W(x) = J_QC/J_N
        double rhoNorm = dy(jo)*dy(jo)/((*stddev_)(jo)*(*stddev_)(jo));
        if (dy(jo) < cleft_) {
          double rhoHuber = (2.0*std::abs(cleft_*dy(jo))-cleft_*cleft_)
            /((*stddev_)(jo)*(*stddev_)(jo));
          (*wghtsqrt_)(jo) = rhoHuber/rhoNorm;
        } else if (dy(jo) > cright_) {
          double rhoHuber = (2.0*cright_*dy(jo)-cright_*cright_)
            /((*stddev_)(jo)*(*stddev_)(jo));
          (*wghtsqrt_)(jo) = rhoHuber/rhoNorm;
        }
      }

      (*wghtsqrt_)(jo) = std::sqrt((*wghtsqrt_)(jo));
    }
  }

  oops::Log::trace() << classname() << "::setupWeights done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVector * ObsError::multiply(const ObsVector & dy) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  ObsVector * res = new ObsVector(dy);
  *res /= *inverseVariance_;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
  return res;
}

// -----------------------------------------------------------------------------

ObsVector * ObsError::inverseMultiply(const ObsVector & dy) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;

  ObsVector * res = new ObsVector(dy);
  *res *= *inverseVariance_;

  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
  return res;
}

// -----------------------------------------------------------------------------

ObsVector * ObsError::multiplyWghtSqrt(const ObsVector & dy) const {
  oops::Log::trace() << classname() << "::multiplyWghtSqrt starting" << std::endl;

  ObsVector * res = new ObsVector(dy);
  if (lvarqc_) {
    *res *= *wghtsqrt_;
  }

  oops::Log::trace() << classname() << "::multiplyWghtSqrt done" << std::endl;
  return res;
}

// -----------------------------------------------------------------------------

void ObsError::randomize(ObsVector & dy) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  dy.random();
  dy *= *stddev_;

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsError::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "ObsError::print not implemeted yet";

  oops::Log::trace() << classname() << "::print end" << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace quench
