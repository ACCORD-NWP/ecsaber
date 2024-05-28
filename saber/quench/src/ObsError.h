/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "oops/interface/ObsAuxControlBase.h"
#include "oops/interface/ObsErrorBase.h"

#include "src/Traits.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ObsSpace;

// -----------------------------------------------------------------------------
/// ObsError class

class ObsError : public oops::ObsErrorBase<Traits> {
  using ObsAuxControlPtrMap_ =
    typename std::map<std::string, std::unique_ptr<oops::ObsAuxControlBase<Traits>>>;

 public:
  static const std::string classname()
    {return "quench::ObsError";}

/// Constructor, destructor
  ObsError(const ObsSpace &,
           const eckit::Configuration &);
  ~ObsError()
    {}

/// Linearize and reset for inner loop
  void linearize(const ObsVec &)
    {}
  void linearize(const ObsAuxControlPtrMap_ &,
                 const ObsVec &dy)
    {}

/// Setup VarQC weights
  void setupWeights(const ObsVec &);

/// Multiply a Departure by \f$R\f$
  ObsVec * multiply(const ObsVec &) const;

/// Multiply a Departure by \f$R^{-1}\f$
  ObsVec * inverseMultiply(const ObsVec &) const;

/// Multiply a Departure by \f$W^1/2\f$
  ObsVec * multiplyWghtSqrt(const ObsVec &) const;

/// Generate random perturbation
  void randomize(ObsVec &) const;

/// Get mean error for Jo table
  double getRMSE() const
    {return stddev_->rms();}

 private:
  void print(std::ostream &) const;

  mutable bool lvarqc_;
  mutable double cleft_;
  mutable double cright_;
  std::unique_ptr<ObsVec> stddev_;
  std::unique_ptr<ObsVec> inverseVariance_;
  std::unique_ptr<ObsVec> wghtsqrt_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
