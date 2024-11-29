/*
 * (C) Copyright 2021- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODELSTATE4D_H_
#define OOPS_GENERIC_PSEUDOMODELSTATE4D_H_

#include <string>
#include <vector>

#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State4D.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic implementation of the pseudo model initialized with 4D State
/// (steps through time by stepping through states in 4D state)
template <typename MODEL>
class PseudoModelState4D : public util::Printable,
                           private eckit::NonCopyable,
                           private util::ObjectCounter<Model<MODEL> > {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;
  typedef State4D<MODEL>           State4D_;

 public:
  static const std::string classname() {return "oops::PseudoModelState4D";}

  /// Initialize pseudo model with \p state4d - 4D state to loop through
  /// in the model run and \p tstep - time resolution of the model
  PseudoModelState4D(const State4D_ & state4d,
                     const util::Duration & tstep = util::Duration(0));

  /// \brief Run the forecast from state \p xx for \p len time, with \p post postprocessors
  /// Does not need to be implemented in the models
  void forecast(State_ & xx, const ModelAux_ &, const util::Duration & len,
                PostProcessor<State_> & post) const;

 private:
  /// Initialize the model (prior to forecast)
  void initialize(State_ &, const ModelAux_ &) const;
  /// Model time step
  void step(State_ &, const ModelAux_ &) const;
  /// Finalize the model (after the forecast)
  void finalize(State_ &) const;
  void print(std::ostream &) const;

  /// Reference to 4D state that is used in the model
  const State4D_ & state4d_;
  /// Model's time resolution
  util::Duration   tstep_;
  /// Index of the current state
  mutable size_t currentstate_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModelState4D<MODEL>::PseudoModelState4D(const State4D_ & state4d,
                                              const util::Duration & tstep)
  : state4d_(state4d), tstep_(tstep) {
  const std::vector<util::DateTime> validTimes = state4d_.times();
  if (validTimes.size() > 1) tstep_ = validTimes[1] - validTimes[0];
  std::cout << tstep_ << std::endl;
  Log::trace() << "PseudoModelState4D<MODEL>::PseudoModelState4D done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void PseudoModelState4D<MODEL>::forecast(State_ &xx, const ModelAux_ &maux,
                                         const util::Duration &len,
                                         PostProcessor<State_> &post) const {
  Log::trace() << "PseudoModelState4D<MODEL>::forecast starting" << std::endl;
  const util::DateTime end(xx.validTime() + len);
  Log::info() << "PseudoModelState4D::forecast: forecast starting: " << xx << std::endl;
  this->initialize(xx, maux);
  post.initialize(xx, end, tstep_);
  size_t step = 0;
  while (xx.validTime() < end) {
    this->step(xx, maux);
    post.process(xx);
    ++step;
  }
  post.finalize(xx);
  this->finalize(xx);
  Log::info() << "PseudoModelState4D::forecast: forecast finished: " << xx << std::endl;
  ASSERT(xx.validTime() == end);

  Log::trace() << "PseudoModelState4D<MODEL>::forecast done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void PseudoModelState4D<MODEL>::initialize(State_ &xx, const ModelAux_ &maux) const {
  Log::trace() << "PseudoModelState4D<MODEL>::initialize starting" << std::endl;
  currentstate_ = 0;
  xx = state4d_[currentstate_];
  Log::trace() << "PseudoModelState4D<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void PseudoModelState4D<MODEL>::step(State_ &xx, const ModelAux_ &maux) const {
  Log::trace() << "PseudoModelState4D<MODEL>::step starting" << std::endl;
  xx = state4d_[currentstate_++];
  Log::trace() << "PseudoModelState4D<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void PseudoModelState4D<MODEL>::finalize(State_ &xx) const {
  Log::trace() << "PseudoModelState4D<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModelState4D<MODEL>::print(std::ostream & os) const {
  os << "Pseudo model stepping through 4D state with " << tstep_ << " time resolution";
}

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODELSTATE4D_H_
