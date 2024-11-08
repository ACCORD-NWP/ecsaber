/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSLOCALIZATIONBASE_H_
#define OOPS_BASE_OBSLOCALIZATIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Printable.h"

namespace oops {

/// Base class for observation-space localization.
/// Defines the interfaces for observation space localization.
/// Use this class as a base class for MODEL-specific implementations.
template<typename MODEL>
class ObsLocalizationBase : public util::Printable,
                                   private boost::noncopyable {
  typedef typename MODEL::GeometryIterator          GeometryIterator_;
  typedef typename MODEL::ObsVector                 ObsVector_;
 public:
  ObsLocalizationBase() = default;
  virtual ~ObsLocalizationBase() = default;

  /// compute obs-space localization: update \p locfactor with observation-space
  /// localization values between observations and \p point in model-space.
  /// update means that locfactor values that are passed in are multiplied by the
  /// locfactor computed inside of this function
  /// Set \p locfactor to missing value for observations that are not local.
  /// Method used in oops. Calls `computeLocalization` abstract method, and
  /// passes MODEL-specific classe to the MODEL-specific
  /// implementation of ObsLocalization.
  void computeLocalization(const GeometryIterator<MODEL> & point,
                           ObsVector<MODEL> & locfactor) const {
    computeLocalization(point.geometryiter(), locfactor.obsvector());
  }

  /// compute obs-space localization: update \p locfactor with observation-space
  /// localization values between observations and \p point in model-space.
  /// Set \p locfactor to missing value for observations that are not local.
  virtual void computeLocalization(const GeometryIterator_ & point,
                                   ObsVector_ & locfactor) const = 0;
};

// =============================================================================

/// ObsLocalization Factory
template <typename MODEL>
class ObsLocalizationFactory {
  typedef ObsSpace<MODEL>  ObsSpace_;
 public:
  static std::unique_ptr<ObsLocalizationBase<MODEL>> create(const eckit::Configuration &,
                                                            const ObsSpace_ &);
 protected:
  explicit ObsLocalizationFactory(const std::string &);
  virtual ~ObsLocalizationFactory() = default;
 private:
  virtual ObsLocalizationBase<MODEL> * make(const eckit::Configuration &,
                                                 const ObsSpace_ &) = 0;
  static std::map < std::string, ObsLocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsLocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsLocalizationMaker : public ObsLocalizationFactory<MODEL> {
  typedef ObsSpace<MODEL>  ObsSpace_;
  virtual ObsLocalizationBase<MODEL> * make(const eckit::Configuration & conf,
                                            const ObsSpace_ & obspace)
    { return new T(conf, obspace.obsspace()); }
 public:
  explicit ObsLocalizationMaker(const std::string & name) :
    ObsLocalizationFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsLocalizationFactory<MODEL>::ObsLocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs localization factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<ObsLocalizationBase<MODEL>> ObsLocalizationFactory<MODEL>::create(
                              const eckit::Configuration & conf, const ObsSpace_ & obspace) {
  Log::trace() << "ObsLocalizationBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization method");
  typename std::map<std::string, ObsLocalizationFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs localization factory.");
  }
  std::unique_ptr<ObsLocalizationBase<MODEL>> ptr(jloc->second->make(conf, obspace));
  Log::trace() << "ObsLocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONBASE_H_
