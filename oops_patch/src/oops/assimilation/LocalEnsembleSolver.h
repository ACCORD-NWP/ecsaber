/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_
#define OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_

#include <Eigen/Dense>
#include <cfloat>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/LocalEnsembleSolverParameters.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/interface/Geometry.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/interface/LinearModel.h"
#include "oops/interface/Model.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObservationSpaces.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observer.h"
#include "oops/interface/State.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
//#include "oops/generic/PseudoLinearModelIncrement4D.h"
#include "oops/generic/PseudoModelState4D.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ModelAuxControl.h"
//#include "oops/interface/ObsDataVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {
  class JediVariables;

/// \brief Base class for LETKF-type solvers
template <typename MODEL>
class LocalEnsembleSolver {
  typedef Departures<MODEL>           Departures_;
  typedef DeparturesEnsemble<MODEL>   DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef typename MODEL::GeometryIterator  GeometryIterator__;
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef ObsAuxControls<MODEL>       ObsAux_;
  typedef ObsAuxIncrements<MODEL>     ObsAuxInc_;
//  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<MODEL>          ObsEnsemble_;
  typedef ObsErrors<MODEL>            ObsErrors_;
  typedef Observations<MODEL>         Observations_;
  typedef ObsLocalizations<MODEL>     ObsLocalizations_;
  typedef ObservationSpaces<MODEL>    ObsSpaces_;
  typedef ObsOperators<MODEL>         ObsOperators_;
  typedef ObsAuxControls<MODEL>       ObsAuxCtrls_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
//  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  typedef State<MODEL>                State_;
  typedef State4D<MODEL>              State4D_;
  typedef Increment<MODEL>            Increment_;
  typedef Increment4D<MODEL>          Increment4D_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;

 public:
  static const std::string classname() {return "oops::LocalEnsembleSolver";}

  /// initialize solver with \p obspaces, \p geometry, full \p config and \p nens ensemble size
  /// \p xbmean state is used if an implementation needs a reference state
  /// solver will use a list of analysis variables specified in \p incvars
  LocalEnsembleSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                      const eckit::Configuration & config, size_t nens, const State4D_ & xbmean,
                      const JediVariables & incvars);
  virtual ~LocalEnsembleSolver() = default;

  /// computes ensemble H(\p xx), returns mean H(\p xx), saves as hofx \p iteration
  virtual Observations_ computeHofX(const StateEnsemble4D_ & xx, size_t iteration,
                      bool readFromDisk, const Model_ & model, const Observations_ & yobs,
                      const ObsAuxCtrls_ & ybias);
  Observations_ computeHofXLinear(const StateEnsemble4D_ & xx, size_t iteration,
                      bool readFromDisk, const Model_ & model, const Observations_ & yobs,
                      const ObsAuxCtrls_ & ybias);
  Observations_ computeHofXNonLinear(const StateEnsemble4D_ & xx, size_t iteration,
                      bool readFromDisk, const Model_ & model, const Observations_ & yobs,
                      const ObsAuxCtrls_ & ybias);

  /// update background ensemble \p bg to analysis ensemble \p for all points on this PE
  virtual void measurementUpdate(const IncrementEnsemble4D_ & bg, IncrementEnsemble4D_ & an);

  /// update background ensemble \p bg to analysis ensemble \p an at a grid point location \p i
  virtual void measurementUpdate(const IncrementEnsemble4D_ & bg,
                                 const GeometryIterator__ & i, IncrementEnsemble4D_ & an) = 0;

  /// copy \p an[\p i] = \p bg[\p i] (e.g. when there are no local observations to update state)
  virtual void copyLocalIncrement(const IncrementEnsemble4D_ & bg,
                                  const GeometryIterator__ & i, IncrementEnsemble4D_ & an) const;

  /// apply posterior inflation to a local ensemble
  void posteriorInflation(const Eigen::MatrixXd & Xb, Eigen::MatrixXd & Xa) const;

  /// compute H(x) based on 4D state \p xx and put the result into \p yy. Also sets up
  /// R_ based on the QC filters run during H(x)
  void computeHofX4DLinear(const eckit::Configuration &, const Model_ &, const StateEnsemble4D_ &,
                           Observations_ &, ObsEnsemble_ &);
  void computeHofX4DNonLinear(const eckit::Configuration &, const Model_ &, const State4D_ &,
                              const ObsAuxCtrls_ &, Observations_ &);
  /// accessor to obs localizations
  const ObsLocalizations_ & obsloc() const {return obsloc_;}
  bool useLinearObserver() { return useLinearObserver_; }

 protected:
  const Geometry_  & geometry_;   ///< Geometry associated with the updated states
  const ObsSpaces_ & obspaces_;   ///< ObsSpaces used in the update
  Departures_ omb_;               ///< obs - mean(H(x)); set in computeHofX method
  DeparturesEnsemble_ Yb_;        ///< ensemble perturbations in the observation space;
                                  ///  set in computeHofX method
  std::unique_ptr<ObsErrors_>   R_;        ///< observation errors, set in computeHofX method
  std::unique_ptr<Departures_> invVarR_;   ///< inverse observation error variance; set in
                                           ///  computeHofX method
  LocalEnsembleSolverParameters options_;

  const State4D_ & xbmean_;     ///< ensemble mean or a control member that will be used to
                                ///  center the prior ensemble
  const JediVariables incvars_;
  const eckit::LocalConfiguration obsconf_;  // configuration for observations
  const eckit::LocalConfiguration observersconf_;  // configuration for observations.observers

 private:
  bool useLinearObserver_;
  ObsLocalizations_ obsloc_;          ///< observation space localization
//  std::vector<ObsDataInt_> qcflags_;  ///< quality control flags
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LocalEnsembleSolver<MODEL>::LocalEnsembleSolver(ObsSpaces_ & obspaces,
                                        const Geometry_ & geometry,
                                        const eckit::Configuration & config, size_t nens,
                                        const State4D_ & xbmean, const JediVariables & incvars)
  : geometry_(geometry),
    obspaces_(obspaces),
    omb_(obspaces_),
    Yb_(obspaces_, nens),
    xbmean_(xbmean),
    incvars_(incvars),
    obsconf_(config, "observations"),
    observersconf_(obsconf_, "ObsTypes"),
    obsloc_(observersconf_, obspaces_) {
  // initialize and print options

  options_.deserialize(config);
  useLinearObserver_ = this->options_.useLinearObserver;
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  Log::info() << "Multiplicative inflation will be applied with multCoeff=" <<
                 inflopt.mult << std::endl;
  if (inflopt.doRtpp()) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                    inflopt.rtpp << std::endl;
  } else {
      Log::info() << "RTPP inflation is not applied rtppCoeff is out of bounds (0,1], rtppCoeff="
                  << inflopt.rtpp << std::endl;
  }
  if (inflopt.doRtps()) {
    Log::info() << "RTPS inflation will be applied with rtpsCoeff=" <<
                    inflopt.rtps << std::endl;
  } else {
    Log::info() << "RTPS inflation is not applied rtpsCoeff is out of bounds (0,1], rtpsCoeff="
                << inflopt.rtps << std::endl;
  }
/*
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObsVector_ qcflags(obspaces_[jj], obspaces_[jj].obsvariables());
    qcflags_.push_back(qcflags);
  }
*/
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalEnsembleSolver<MODEL>::measurementUpdate
        (const IncrementEnsemble4D_ & bg, IncrementEnsemble4D_ & an) {
    for (auto i = geometry_.geometry().begin(); i != geometry_.geometry().end(); ++i) {
      measurementUpdate(bg, i, an);
    }
}
// -----------------------------------------------------------------------------

template <typename MODEL>
Observations<MODEL> LocalEnsembleSolver<MODEL>::computeHofX(const StateEnsemble4D_ & ens_xx,
  size_t iteration, bool readFromDisk, const Model_ & model, const Observations_ & yobs,
  const ObsAuxCtrls_ & ybias) {
  util::Timer timer(classname(), "computeHofX");

  Observations_ yb_mean(obspaces_);

  if (useLinearObserver()) {
    yb_mean = computeHofXLinear(ens_xx, iteration, readFromDisk, model, yobs, ybias);
  } else {
    yb_mean = computeHofXNonLinear(ens_xx, iteration, readFromDisk, model, yobs, ybias);
  }

  // return mean H(x)
  return yb_mean;
}
// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalEnsembleSolver<MODEL>::computeHofX4DLinear(const eckit::Configuration & config,
                                                     const Model_ & model,
                                                     const StateEnsemble4D_ & xx,
                                                     Observations_ & yy_mean,
                                                     ObsEnsemble_ & yy) {
/*
  ModelAux_ moderr(geometry_, model, eckit::LocalConfiguration());
  ModelAuxInc_  moderrinc(geometry_, eckit::LocalConfiguration());
  ObsAux_  obsaux(obspaces_, obsconf_);
  ObsAuxInc_  obsauxinc(obspaces_, obsconf_);

  // compute forecast length from State4D times
  const std::vector<util::DateTime> times = xx[0].times();
  const util::Duration flength = times[times.size()-1] - times[0];
  // default_tstep = 2*observation window is passed to PseudoModel as the default
  // pseudomodel time step. It is only used when State4D has a single state, to enable
  // processing of all observations in the specified window regardless of where in
  // the time window the state is. Observations in
  // ( max(winbgn, xx.time - tstep/2); min(winend, xx.time + tstep/2) ] are
  // processed in H(x).
  const util::Duration default_tstep = (obspaces_.windowEnd() - obspaces_.windowStart()) * 2;
  R_.reset(new ObsErrors_(obspaces_));

  // Setup pseudo model to run on ensemble mean
  State_ init_xx = xbmean_[0];
  std::unique_ptr<PseudoModel_> pseudomodel(new PseudoModel_(xbmean_, default_tstep));
  const Model_ model(std::move(pseudomodel));

  // setup postprocessors and nonlinear observers for the "nonlinear" model run on the mean
  PostProcessor<State_> post;
  PostProcessorTLAD<MODEL> posttraj;
  Observers_ hofx(obspaces_, obsconf_);

  // setup postprocessors and linear observers for the "linear" model run on the ensemble
  // perturbations
  PostProcessor<Increment_> posttl;
  PostProcessorTLAD<MODEL> posttrajtl;
  ObserversTLAD_ linear_hofx(obspaces_, obsconf_);

  // initialize nonlinear model postprocessor
  hofx.initialize(geometry_, obsaux, *R_, post, config);

  // add linearized H(x) to the nonlinear model postprocessor
  linear_hofx.initializeTraj(geometry_, obsaux, posttraj);
  // create TrajectorySaver with hofx_linear, and enroll in post
  post.enrollProcessor(new TrajectorySaver<MODEL>(eckit::LocalConfiguration(),
                                                  geometry_, posttraj));

  // run nonlinear model on the ensemble mean
  model.forecast(init_xx, moderr, flength, post);

  // compute nonlinear H(x_mean)
  hofx.finalize(yy_mean, qcflags_);
  linear_hofx.finalizeTraj(qcflags_);

  // add linearized H(x) to the linear model postprocessor
  linear_hofx.initializeTL(posttrajtl);
  for (size_t jens = 0; jens < xx.size(); ++jens) {
    // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
    Increment4D_ dx(geometry_, xx[jens].variables(), times);
    dx.diff(xx[jens], xbmean_);
    Increment_ init_dx = dx[0];
    std::unique_ptr<PseudoLinearModel_> pseudolinearmodel =
         std::make_unique<PseudoLinearModel_>(dx, default_tstep);
    const LinearModel_ linear_model(std::move(pseudolinearmodel));
    // run linear model on the ensemble perturbation, compute linear H*dx
    linear_model.forecastTL(init_dx, moderrinc, flength, posttl, posttrajtl);
    linear_hofx.finalizeTL(obsauxinc, Yb_[jens]);
    yy[jens] = yy_mean;
    yy[jens] += Yb_[jens];
  }
*/
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Observations<MODEL> LocalEnsembleSolver<MODEL>::computeHofXLinear(
                                                   const StateEnsemble4D_ & ens_xx,
                                                   size_t iteration,
                                                   bool readFromDisk,
                                                   const Model_ & model,
                                                   const Observations_ & yobs,
                                                   const ObsAuxCtrls_ & ybias) {
  util::Timer timer(classname(), "computeHofXLinear");

  ASSERT(ens_xx.size() == Yb_.size());

  const size_t nens = ens_xx.size();
  ObsEnsemble_ obsens(obspaces_, nens);
  Observations_ y_mean_xb(obspaces_);

  if (readFromDisk) {
    // read hofx from disk
    for (size_t jj = 0; jj < nens; ++jj) {
      eckit::LocalConfiguration conf;
      conf.set("ObsValue", "hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
      obsens[jj].read(conf);
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
    }
    R_.reset(new ObsErrors_(obspaces_));
    eckit::LocalConfiguration conf;
    conf.set("ObsValue", "hofx_y_mean_xb"+std::to_string(iteration));
    y_mean_xb.read(conf);
  } else {
    // compute and save H(x)

    // save QC filters and ob errors to be used for all other members
    // do not save H(X) (saved explicitly below)
    eckit::LocalConfiguration config;

    // save hofx means that hofx will be written out into ObsSpace;
    // if run computeHofX4D several times with save hofx on,
    // the hofx will be overwritten,
    // unless each time specifying iteration differently in the passed config.
    // TODO(someone): this needs to be revisited so the flags are saved for the mean.
    config.set("save hofx", false);
    config.set("save qc", false);
    config.set("save obs errors", false);
    config.set("iteration", std::to_string(iteration));

    computeHofX4DLinear(config, model, ens_xx, y_mean_xb, obsens);
    for (size_t jj = 0; jj < nens; ++jj) {
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
      eckit::LocalConfiguration conf;
      conf.set("ObsValue", "hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
      obsens[jj].save(conf);
    }

    // Compute H(mean(Xb))
    // set QC for the mean
    config.set("save qc", true);
    config.set("save obs errors", true);

    eckit::LocalConfiguration conf;
    conf.set("ObsValue", "hofx_y_mean_xb"+std::to_string(iteration));
    y_mean_xb.save(conf);

    // QC flags and Obs errors are set to that of the H(mean(Xb))
    R_->save("ObsError");
  }

  // set inverse variances
  invVarR_.reset(new Departures_(obspaces_));
  R_->inverseVariance(*invVarR_);

  // calculate H(x) ensemble mean
  Observations_ yb_mean(obsens.mean());

  // treat the special case of nens=1
  // default option: xbmean_=mean(xb) then yb_mean == y_mean_xb and action below is a tautology
  // if use control member==true: xbmean_ was read from the controll member,
  //                              then using H(xbmean_) is expected by downstream applications
  if (nens == 1) {yb_mean = y_mean_xb;}

  // mask H(x) ensemble perturbations
  for (size_t iens = 0; iens < nens; ++iens) {
    if (readFromDisk) {
      for (size_t jj = 0; jj < Yb_[iens].size(); ++jj) {
        Yb_[iens][jj] = obsens[iens][jj];
        Yb_[iens][jj] -= yb_mean[jj];
      }
    }
    invVarR_->mask(Yb_[iens]);
    Yb_[iens].mask(*invVarR_);
  }

  // calculate obs departures and mask with qc flag
  for (size_t jj = 0; jj < omb_.size(); ++jj) {
    omb_[jj] = yobs[jj];
    omb_[jj] -= yb_mean[jj];
  }
  invVarR_->mask(omb_);
  omb_.mask(*invVarR_);

  // return mean H(x)
  return yb_mean;
}
// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalEnsembleSolver<MODEL>::computeHofX4DNonLinear(const eckit::Configuration & config,
                                                        const Model_ & model,
                                                        const State4D_ & xx,
                                                        const ObsAuxCtrls_ & ybias,
                                                        Observations_ & yy) {
  // compute forecast length from State4D times
  const std::vector<util::DateTime> times = xx.times();
  const util::Duration flength = times[times.size()-1] - times[0];
  // default_tstep = 2*observation window is passed to PseudoModel as the default
  // pseudomodel time step. It is only used when State4D has a single state, to enable
  // processing of all observations in the specified window regardless of where in
  // the time window the state is. Observations in
  // ( max(winbgn, xx.time - tstep/2); min(winend, xx.time + tstep/2) ] are
  // processed in H(x).
  const util::Duration default_tstep = (obspaces_.windowEnd() - obspaces_.windowStart()) * 2;
  // Setup PseudoModelState4D
  PseudoModel_ pseudomodel(xx, default_tstep);
  // Setup and run the model forecast with observers
  State_ init_xx = xx[0];
  PostProcessor<State_> post;
  ObsOperators_ hop(obspaces_);
  ModelAux_ moderr(geometry_, model, eckit::LocalConfiguration());
  ModelAuxInc_  moderrinc(geometry_, eckit::LocalConfiguration());
  std::shared_ptr<Observer<MODEL, State_> > pobs(
    new Observer<MODEL, State_>(obspaces_, hop, ybias));
  post.enrollProcessor(pobs);
  pseudomodel.forecast(init_xx, moderr, flength, post);
  yy = *pobs->release();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Observations<MODEL> LocalEnsembleSolver<MODEL>::computeHofXNonLinear(
                                                   const StateEnsemble4D_ & ens_xx,
                                                   size_t iteration,
                                                   bool readFromDisk,
                                                   const Model_ & model,
                                                   const Observations_ & yobs,
                                                   const ObsAuxCtrls_ & ybias) {
  util::Timer timer(classname(), "computeHofXNonLinear");

  ASSERT(ens_xx.size() == Yb_.size());

  const size_t nens = ens_xx.size();
  ObsEnsemble_ obsens(obspaces_, nens);
  Observations_ y_mean_xb(obspaces_);

  if (readFromDisk) {
    // read hofx from disk
    for (size_t jj = 0; jj < nens; ++jj) {
      eckit::LocalConfiguration conf;
      conf.set("ObsValue", "hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
      obsens[jj].read(conf);
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
    }
    R_.reset(new ObsErrors_(obspaces_));
    eckit::LocalConfiguration conf;
    conf.set("ObsValue", "hofx_y_mean_xb"+std::to_string(iteration));
    y_mean_xb.read(conf);
  } else {
    // compute and save H(x)

    // save QC filters and ob errors to be used for all other members
    // do not save H(X) (saved explicitly below)
    eckit::LocalConfiguration config;

    // save hofx means that hofx will be written out into ObsSpace;
    // if run computeHofX4D several times with save hofx on,
    // the hofx will be overwritten,
    // unless each time specifying iteration differently in the passed config.
    config.set("save hofx", false);
    config.set("save qc", false);
    config.set("save obs errors", false);
    config.set("iteration", std::to_string(iteration));

    for (size_t jj = 0; jj < nens; ++jj) {
      computeHofX4DNonLinear(config, model, ens_xx[jj], ybias, obsens[jj]);
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
//      obsens[jj].save("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
    }

    // Compute H(mean(Xb))
    // set QC for the mean
    config.set("save qc", true);
    config.set("save obs errors", true);

    computeHofX4DNonLinear(config, model, xbmean_, ybias, y_mean_xb);

    // Setup model and obs biases; obs errors
    R_.reset(new ObsErrors_(obspaces_));
    R_->linearize(yobs);
//    y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));

    // QC flags and Obs errors are set to that of the H(mean(Xb))
//    R_->save("ObsError");
  }

  // set inverse variances
  invVarR_.reset(new Departures_(obspaces_));
  R_->inverseVariance(*invVarR_);

  // calculate H(x) ensemble mean
  Observations_ yb_mean(obsens.mean());

  // treat the special case of nens=1
  // default option: xbmean_=mean(xb) then yb_mean == y_mean_xb and action below is a tautology
  // if use control member==true: xbmean_ was read from the controll member,
  //                              then using H(xbmean_) is expected by downstream applications
  if (nens == 1) {yb_mean = y_mean_xb;}

  // calculate H(x) ensemble perturbations
  for (size_t iens = 0; iens < nens; ++iens) {
    if (readFromDisk) {
      for (size_t jj = 0; jj < Yb_[iens].size(); ++jj) {
        Yb_[iens][jj] = obsens[iens][jj];
        Yb_[iens][jj] -= yb_mean[jj];
      }
    }
    invVarR_->mask(Yb_[iens]);
    Yb_[iens].mask(*invVarR_);
  }

  // calculate obs departures and mask with qc flag
  for (size_t jj = 0; jj < omb_.size(); ++jj) {
    omb_[jj] = yobs[jj];
    omb_[jj] -= yb_mean[jj];
  }
  invVarR_->mask(omb_);
  omb_.mask(*invVarR_);

  // return mean H(x)
  return yb_mean;
}


// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalEnsembleSolver<MODEL>::copyLocalIncrement(const IncrementEnsemble4D_ & bkg_pert,
                                                    const GeometryIterator__ & i,
                                                    IncrementEnsemble4D_ & ana_pert) const {
  // ana_pert[i]=bkg_pert[i]
  for (size_t itime=bkg_pert[0].first(); itime < bkg_pert[0].last()+1; ++itime) {
    for (size_t iens=0; iens < bkg_pert.size(); ++iens) {
      LocalIncrement gp = bkg_pert[iens][itime].increment().getLocal(i);
      ana_pert[iens][itime].increment().setLocal(gp, i);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalEnsembleSolver<MODEL>::posteriorInflation(
                                  const Eigen::MatrixXd & Xb, Eigen::MatrixXd & Xa) const {
    const size_t nens = Xa.cols();
    const LocalEnsembleSolverInflationParameters & inflopt = options_.infl;

    // RTPP inflation
    if (inflopt.doRtpp()) {
      Xa = (1-inflopt.rtpp)*Xa+inflopt.rtpp*Xb;
    }

    // RTPS inflation
    const double eps = DBL_EPSILON;
    if (inflopt.doRtps()) {
      // posterior spread
      Eigen::ArrayXd asprd = Xa.array().square().rowwise().sum()/(nens-1);
      asprd = asprd.sqrt();
      asprd = (asprd < eps).select(eps, asprd);  // avoid nan overflow for vars with no spread

      // prior spread
      Eigen::ArrayXd fsprd = Xb.array().square().rowwise().sum()/(nens-1);
      fsprd = fsprd.sqrt();
      fsprd = (fsprd < eps).select(eps, fsprd);

      // rtps inflation factor
      Eigen::ArrayXd rtpsInfl = inflopt.rtps*((fsprd-asprd)/asprd) + 1;
      rtpsInfl = (rtpsInfl < inflopt.rtpsInflMin()).select(inflopt.rtpsInflMin(), rtpsInfl);
      rtpsInfl = (rtpsInfl > inflopt.rtpsInflMax()).select(inflopt.rtpsInflMax(), rtpsInfl);

      // inflate perturbation matrix
      Xa.array().colwise() *= rtpsInfl;
    }
}

// =============================================================================

/// \brief factory for LETKF solvers
template <typename MODEL>
class LocalEnsembleSolverFactory {
  typedef Geometry<MODEL>           Geometry_;
  typedef ObservationSpaces<MODEL>  ObsSpaces_;
  typedef State4D<MODEL>            State4D_;
 public:
  static std::unique_ptr<LocalEnsembleSolver<MODEL>> create(ObsSpaces_ &, const Geometry_ &,
                                                        const eckit::Configuration &,
                                                        size_t, const State4D_ &,
                                                        const JediVariables &);
  virtual ~LocalEnsembleSolverFactory() = default;
 protected:
  explicit LocalEnsembleSolverFactory(const std::string &);
 private:
  virtual LocalEnsembleSolver<MODEL> * make(ObsSpaces_ &, const Geometry_ &,
                                        const eckit::Configuration &, size_t,
                                        const State4D_ &, const JediVariables &) = 0;
  static std::map < std::string, LocalEnsembleSolverFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalEnsembleSolverFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LocalEnsembleSolverMaker : public LocalEnsembleSolverFactory<MODEL> {
  typedef Geometry<MODEL>           Geometry_;
  typedef ObservationSpaces<MODEL>  ObsSpaces_;
  typedef State4D<MODEL>            State4D_;

  virtual LocalEnsembleSolver<MODEL> * make(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                        const eckit::Configuration & conf, size_t nens,
                                        const State4D_ & xbmean, const JediVariables & incvars)
    { return new T(obspaces, geometry, conf, nens, xbmean, incvars); }
 public:
  explicit LocalEnsembleSolverMaker(const std::string & name)
    : LocalEnsembleSolverFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
LocalEnsembleSolverFactory<MODEL>::LocalEnsembleSolverFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in local ensemble solver factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LocalEnsembleSolver<MODEL>>
LocalEnsembleSolverFactory<MODEL>::create(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                  const eckit::Configuration & conf, size_t nens,
                                  const State4D_ & xbmean, const JediVariables & incvars) {
  Log::trace() << "LocalEnsembleSolver<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("local ensemble DA.solver");
  typename std::map<std::string, LocalEnsembleSolverFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in local ensemble solver factory." << std::endl;
    Log::error() << "LETKF solver Factory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LocalEnsembleSolverFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " LocalEnsembleSolver" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in local ensemble solver factory.");
  }
  std::unique_ptr<LocalEnsembleSolver<MODEL>>
    ptr(jloc->second->make(obspaces, geometry, conf, nens, xbmean, incvars));
  Log::trace() << "LocalEnsembleSolver<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_
