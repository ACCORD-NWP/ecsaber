/*
 * (C) Copyright 2021-2023 UCAR
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/IncrCtlVec.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

#include "util/abor1_cpp.h"
#include "util/Logger.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance4D : public oops::ModelSpaceCovariance4DBase<MODEL> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef oops::Increment4D<MODEL>                             Increment4D_;
  typedef oops::IncrCtlVec<MODEL>                              IncrCtlVec_;
  typedef oops::IncrEnsCtlVec<MODEL>                           IncrEnsCtlVec_;
  typedef oops::IncrModCtlVec<MODEL>                           IncrModCtlVec_;
  typedef oops::LinearVariableChange<MODEL>                    LinearVariableChange_;
  typedef oops::State4D<MODEL>                                 State4D_;
  typedef oops::Variables<MODEL>                               Variables_;

 public:
  static const std::string classname() {return "saber::ErrorCovariance4D";}

  ErrorCovariance4D(const Geometry_ &, const Variables_ &, const eckit::Configuration &,
                    const State4D_ &);
  ~ErrorCovariance4D();

  // Methods
  void advectedLinearize(const State4D_ &, const Geometry_ &,
                         const eckit::Configuration &) override;
  void advectedMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void advectedInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void advectedMultiplySqrt(const IncrCtlVec_ &, Increment4D_ &) const override;
  void advectedMultiplySqrtTrans(const Increment4D_ &, IncrCtlVec_ &) const override;
  void randomize(Increment4D_ &) const override;
  const oops::ModelSpaceCovarianceBase<MODEL> &covar() const { return *cov3d_; }

  // Control Vector
  IncrModCtlVec_ *newIncrModCtlVec() const override {
    return new IncrModCtlVec_(this->ctlVecSize());
  }
  IncrEnsCtlVec_ *newIncrEnsCtlVec() const override {
    return new IncrEnsCtlVec_();
  }

  size_t ctlVecSize() const
    {return blockChain_->ctlVecSize();}   // TODO(Benjamin): necessary?

 private:
  void print(std::ostream &) const;

  /// Chain of blocks (hybrid or ensemble or parametric)
  std::unique_ptr<SaberBlockChainBase> blockChain_;

  /// Geometry UID
  std::string uid_ = "";

  /// Linear variable change
  std::unique_ptr<Variables_> BVars_;
  std::unique_ptr<Variables_> anaVars_;
  std::vector<std::unique_ptr<LinearVariableChange_>> linVarChg_;

  // Dummy 3D covariance
  std::unique_ptr<oops::ModelSpaceCovarianceBase<MODEL>> cov3d_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance4D<MODEL>::ErrorCovariance4D(const Geometry_ & geom,
                                            const Variables_ & incVars,
                                            const eckit::Configuration & config,
                                            const State4D_ & xb)
  : oops::ModelSpaceCovariance4DBase<MODEL>::ModelSpaceCovariance4DBase(geom, config)
{
  oops::Log::trace() << "ErrorCovariance4D::ErrorCovariance4D starting" << std::endl;

  // Deserialize parameters
  ErrorCovarianceParameters params;
  params.deserialize(config);

  const eckit::LocalConfiguration cvconf = config.getSubConfiguration("linear variable change");
  if (!cvconf.empty()) {
     // Setup linear variable changes
    if (cvconf.has("input variables")) {
      eckit::LocalConfiguration inputVars;
      inputVars.set("variables", cvconf.getStringVector("input variables"));
      BVars_.reset(new Variables_(inputVars));
    }
    eckit::LocalConfiguration outputVars;
    if (cvconf.has("output variables")) {
      outputVars.set("variables", cvconf.getStringVector("output variables"));
    }
    anaVars_.reset(new Variables_(outputVars));
    for (size_t jt = 0; jt < xb.statesNumber(); ++jt) {
      linVarChg_.push_back(std::make_unique<LinearVariableChange_>(geom, cvconf));
      linVarChg_[jt]->changeVarTraj(xb[jt], *anaVars_);
    }
  } else {
    // Setup B and analysis variables
    BVars_.reset(new Variables_(incVars));
    anaVars_.reset(new Variables_(incVars));
  }


  oops::Log::trace() << "ErrorCovariance4D::ErrorCovariance4D done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance4D<MODEL>::~ErrorCovariance4D() {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::~ErrorCovariance4D starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance4D");
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::~ErrorCovariance4D done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance4D<MODEL>::advectedLinearize(const State4D_ & xb,
                                                 const Geometry_ & geom,
                                                 const eckit::Configuration & config) {
  oops::Log::trace() << "ErrorCovariance4D::advectedLinearize starting" << std::endl;

  // Deserialize parameters
  ErrorCovarianceParameters params;
  params.deserialize(config);

  // Check whether a new setup is required
  bool newSetup = true;

  // Condition based on geometry
  const std::string uid = util::getGridUid(geom.geometry().functionSpace());
  if (uid_ != "") {
    if (uid_ == uid) {
      newSetup = false;
    }
  }
  uid_ = uid;

  if (newSetup) {
    // JEDI compatibility
    State4D_ fg(xb);

    // Reset components
    blockChain_.reset();

    // Local copy of background and first guess that can undergo interpolation
    std::unique_ptr<oops::FieldSet4D> fset4dXb;
    std::unique_ptr<oops::FieldSet4D> fset4dFg;

    // Change resolution if needed
    if (params.changeBackgroundResolution) {
      const State4D_ xb_lowres(geom, xb);
      const State4D_ fg_lowres(geom, fg);
      const oops::FieldSet4D fset4dXbTmp(xb_lowres);
      const oops::FieldSet4D fset4dFgTmp(fg_lowres);
      fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
      fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
    } else {
      const oops::FieldSet4D fset4dXbTmp(xb);
      const oops::FieldSet4D fset4dFgTmp(fg);
      fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
      fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
    }

    // Initialize outer variables
    const std::vector<std::size_t> vlevs = geom.geometry().variableSizes(BVars_->variables());
    oops::JediVariables outerVars(BVars_->variables().variablesList());
    for (std::size_t i = 0; i < vlevs.size() ; ++i) {
      outerVars[i].setLevels(vlevs[i]);
    }

    // Create blockchain
    blockChain_ = SaberBlockChainFactory<MODEL>::create(
          geom,
          outerVars,
          *fset4dXb,
          *fset4dFg,
          config);
  }

  oops::Log::trace() << "ErrorCovariance4D advectedLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4D<MODEL>::randomize(Increment4D_ & dx) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // This extra fieldset is only needed for backward compatibility in tests:
  // this allows the fields in the final fieldset to be in the same order as
  // blockChain_->outerVariables()
  oops::FieldSet4D fset4dSum(dx.times(), eckit::mpi::self(), dx.geometry().geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
    fset4dSum[jtime].init(blockChain_->outerFunctionSpace(),
                          blockChain_->outerVariables(),
                          0.0);
  }

  // Create FieldSet4D, run randomize on it
  oops::FieldSet4D fset4d(dx.times(), eckit::mpi::self(), dx.geometry().geometry().getComm());

  // Draw a random sample from the covariance
  blockChain_->randomize(fset4d);

  // For backward compatibility in tests
  fset4dSum += fset4d;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dx.times().size(); ++jtime) {
    dx[jtime].increment().fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  if (linVarChg_.size() > 0) {
    // Apply linear variable changes
    for (size_t jtime = 0; jtime <= dx.times().size(); ++jtime) {
      linVarChg_[jtime]->changeVarTL(dx[jtime], *anaVars_);
    }
  }

  oops::Log::trace() << "ErrorCovariance4D<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance4D<MODEL>::advectedMultiply(const Increment4D_ &dxi,
                                                Increment4D_ &dxo) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiply starting" << std::endl;
  util::Timer timer(classname(), "advectedMultiply");

  dxo = dxi;
  if (linVarChg_.size() > 0) {
    // Apply adjoint variable change (to control variables)
    for (size_t jtime = 0; jtime <= dxi.times().size(); ++jtime) {
      linVarChg_[jtime]->changeVarAD(dxo[jtime], *BVars_);
    }
  }
  oops::FieldSet4D fset4d(dxo);

  // Extra fieldset only needed for backward compatibility in tests
  // that allows the fields in the final fieldset to be in the same order as
  // blockChain_->outerVariables()
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4d);
  fset4dSum.zero();

  // Apply covariance multiplication
  blockChain_->multiply(fset4d);

  // For backward compatibility in tests
  fset4dSum += fset4d;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dxo.times().size(); ++jtime) {
    dxo[jtime].increment().fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  if (linVarChg_.size() > 0) {
    // Apply control to analysis/model variable change
    for (size_t jtime = 0; jtime <= dxo.times().size(); ++jtime) {
      linVarChg_[jtime]->changeVarTL(dxo[jtime], *anaVars_);
    }
  }

  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance4D<MODEL>::advectedInverseMultiply(const Increment4D_ &dxi,
                                                       Increment4D_ &dxo) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "advectedInverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance4D<MODEL>::advectedMultiplySqrt(const IncrCtlVec_ &dv,
                                                    Increment4D_ &dx) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiplySqrt starting" << std::endl;

  // This extra fieldset is only needed for backward compatibility in tests:
  // this allows the fields in the final fieldset to be in the same order as
  // blockChain_->outerVariables()
  oops::FieldSet4D fset4dSum(dx.times(), eckit::mpi::self(), dx.geometry().geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
    fset4dSum[jtime].init(blockChain_->outerFunctionSpace(),
                          blockChain_->outerVariables(),
                          0.0);
  }

  // Create FieldSet4D, run square-root on it
  oops::FieldSet4D fset4d(dx.times(), eckit::mpi::self(), dx.geometry().geometry().getComm());

  // Apply blockchain square-root
  blockChain_->multiplySqrt(dv.modCtlVec().genCtlVec().data(), fset4d, 0);

  // For backward compatibility in tests
  fset4dSum += fset4d;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dx.times().size(); ++jtime) {
    dx[jtime].increment().fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  if (linVarChg_.size() > 0) {
    // Apply linear variable changes
    for (size_t jtime = 0; jtime <= dx.times().size(); ++jtime) {
      linVarChg_[jtime]->changeVarTL(dx[jtime], *anaVars_);
    }
  }

  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance4D<MODEL>::advectedMultiplySqrtTrans(const Increment4D_ &dx,
                                                         IncrCtlVec_ &dv) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiplySqrtTrans starting" << std::endl;

  // Setup FieldSet4D
  oops::FieldSet4D fset4d(dx);
  if (linVarChg_.size() > 0) {
    // Copy input increment and apply adjoint variable change (to control variables)
    Increment4D_ dxTmp(dx);
    for (size_t jtime = 0; jtime <= dx.times().size(); ++jtime) {
      linVarChg_[jtime]->changeVarAD(dxTmp[jtime], *BVars_);
    }
    for (size_t jtime = 0; jtime <= dx.times().size(); ++jtime) {
      fset4d[jtime].deepCopy(dxTmp[jtime].increment().fieldSet());
    }
  }

  // Apply blockchain square-root adjoint
  blockChain_->multiplySqrtAD(fset4d, dv.modCtlVec().genCtlVec().data(), 0);

  oops::Log::trace() << "ErrorCovariance4D<MODEL>::advectedMultiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4D<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance4D<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance4D<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------
// ErrorCovariance as a subcase of ErrorCovariance4D
// -----------------------------------------------------------------------------
template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL> {
  using ErrorCovariance4D_ = ErrorCovariance4D<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using IncrCtlVec_ = oops::IncrCtlVec<MODEL>;
  using IncrEnsCtlVec_ = oops::IncrEnsCtlVec<MODEL>;
  using IncrModCtlVec_ = oops::IncrModCtlVec<MODEL>;
  using State_ = oops::State<MODEL>;
  using State4D_ = oops::State4D<MODEL>;
  using Variables_ = oops::Variables<MODEL>;

 public:
  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const Variables_ &, const eckit::Configuration &,
                  const State_ &);
  ~ErrorCovariance();

  // Methods
  void linearize(const State_ &, const Geometry_ &, const eckit::Configuration &) override;
  void multiply(const Increment_ &, Increment_ &) const override;
  void inverseMultiply(const Increment_ &, Increment_ &) const override;
  void multiplySqrt(const IncrCtlVec_ &, Increment_ &) const override;
  void multiplySqrtTrans(const Increment_ &, IncrCtlVec_ &) const override;
  void randomize(Increment_ &) const override;

  // Control Vector
  IncrModCtlVec_ *newIncrModCtlVec() const override {
    return new IncrModCtlVec_(Bmat4D_->ctlVecSize());
  }
  IncrEnsCtlVec_ *newIncrEnsCtlVec() const override {
    return new IncrEnsCtlVec_();
  }

 private:
  void print(std::ostream &) const;

  /// ErrorCovariance 4D
  std::unique_ptr<ErrorCovariance4D_> Bmat4D_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const Variables_ & incVars,
                                        const eckit::Configuration & config,
                                        const State_ & xb3D)
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // 4D compatibility
  State4D_ xb;
  xb.push_back(xb3D);

  // ErrorCovariance4D setup
  Bmat4D_.reset(new ErrorCovariance4D_(geom, incVars, config, xb));

  oops::Log::trace() << "ErrorCovariance::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::linearize(const State_ & xb3D,
                                       const Geometry_ & geom,
                                       const eckit::Configuration & config) {
  // 4D compatibility
  State4D_ xb;
  xb.push_back(xb3D);

  // ErrorCovariance4D linearize
  Bmat4D_->linearize(xb, geom, config);

  oops::Log::trace() << "ErrorCovariance linearized." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::randomize(Increment_ & dx3d) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // 4D compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});

  // ErrorCovariance4D randomize
  Bmat4D_->randomize(dx);

  // 4D compatibility
  dx3d = dx[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiply(const Increment_ &dx3di, Increment_ &dx3do) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // 4D compatibility
  Increment4D_ dxi(dx3di.geometry(), dx3di.variables(), {dx3di.validTime()});
  Increment4D_ dxo(dx3do.geometry(), dx3do.variables(), {dx3do.validTime()});
  dxi[0] = dx3di;

  // ErrorCovariance4D multiply
  Bmat4D_->advectedMultiply(dxi, dxo);

  // 4D compatibility
  dx3do = dxo[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::inverseMultiply(const Increment_ &dx3di, Increment_ &dx3do) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");

  // 4D compatibility
  Increment4D_ dxi(dx3di.geometry(), dx3di.variables(), {dx3di.validTime()});
  Increment4D_ dxo(dx3do.geometry(), dx3do.variables(), {dx3do.validTime()});
  dxi[0] = dx3di;

  // ErrorCovariance4D inverse multiply
  Bmat4D_->advectedInverseMultiply(dxi, dxo);

  // 4D compatibility
  dx3do = dxo[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiplySqrt(const IncrCtlVec_ &dv,
                                          Increment_ &dx3d) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrt starting" << std::endl;

  // 4D compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});

  // ErrorCovariance4D square-root multiply
  Bmat4D_->advectedMultiplySqrt(dv, dx);

  // 4D compatibility
  dx3d = dx[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiplySqrtTrans(const Increment_ &dx3d,
                                               IncrCtlVec_ &dv) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrtTrans starting" << std::endl;

  // 4D compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});
  dx[0] = dx3d;

  // ErrorCovariance4D transposed square-root multiply
  Bmat4D_->advectedMultiplySqrtTrans(dx, dv);

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
