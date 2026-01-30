/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/FieldSets.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/generic/GenericCtlVec.h"
#include "oops/generic/LocalizationBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberParametricBlockChain.h"
#include "oops/util/ECUtilities.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template<typename MODEL>
class Localization : public oops::LocalizationBase<MODEL> {
  typedef oops::GenericCtlVec                 GenericCtlVec_;
  typedef oops::Geometry<MODEL>               Geometry_;
  typedef oops::Increment<MODEL>              Increment_;
  typedef oops::Model<MODEL>                  Model_;
  typedef oops::Variables<MODEL>              Variables_;

 public:
  Localization(const Geometry_ &,
               const Variables_ &,
               const eckit::Configuration &);
  ~Localization();
  static const std::string classname() {return "saber::Localization";}

  void multiply(Increment_ &) const override;

  // Square-root formulation
  size_t ctlVecSize() const {return loc_->ctlVecSize();}
  void multiplySqrt(const GenericCtlVec_ & dv, Increment_ & dx) const;
  void multiplySqrtTrans(const Increment_ & dx, GenericCtlVec_ & dv) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SaberParametricBlockChain> loc_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geom,
                                  const Variables_ & incVarsNoMeta,
                                  const eckit::Configuration & conf)
  : loc_()
{
  oops::Log::trace() << "Localization::Localization starting" << std::endl;
  util::Timer timer(classname(), "Localization");

  // Create dummy time
  util::DateTime dummyTime(1977, 5, 25, 0, 0, 0);

  // Initialize increment variables with levels metadata
  const std::vector<std::size_t> vlevs = geom.geometry().variableSizes(incVarsNoMeta.variables());
  oops::JediVariables incVars(incVarsNoMeta.variables().variablesList());
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    incVars[i].setLevels(vlevs[i]);
  }

  // Create dummy xb and fg
  oops::FieldSet3D fset3d(dummyTime, eckit::mpi::comm());
  fset3d.deepCopy(util::createFieldSet(geom.geometry().functionSpace(), incVars, 0.0));
  oops::FieldSet4D xb4d(fset3d);
  oops::FieldSet4D fg4d(fset3d);

  // 3D localization always used here (4D aspects handled in oops::Localization),
  // so this parameter can be anything.
  eckit::LocalConfiguration confUpdated(conf);
  confUpdated.set("time covariance", "univariate");

  // Initialize localization blockchain
  loc_ = std::make_unique<SaberParametricBlockChain>(geom,
              incVars, xb4d, fg4d, conf);

  oops::Log::trace() << "Localization:Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  util::Timer timer(classname(), "~Localization");
  oops::Log::trace() << "Localization:~Localization destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiplySqrt(const GenericCtlVec_ & dv,
                                       Increment_ & dx) const {
  util::Timer timer(classname(), "multiplySqrt");
  oops::Log::trace() << "Localization:multiplySqrt starting" << std::endl;

  // SABER block chain square-root
  oops::FieldSet3D fset3d(dx.validTime(), eckit::mpi::comm());
  oops::FieldSet4D fset4d(fset3d);
  loc_->multiplySqrt(dv.data(), fset4d, 0);

  // ATLAS fieldset to Increment_
  dx.increment().fromFieldSet(fset4d[0].fieldSet());

  oops::Log::trace() << "Localization:multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiplySqrtTrans(const Increment_ & dx,
                                            GenericCtlVec_ & dv) const {
  util::Timer timer(classname(), "multiplySqrtTrans");
  oops::Log::trace() << "Localization:multiplySqrtTrans starting" << std::endl;

  // SABER block chain square-root adjoint
  oops::FieldSet3D fset3d(dx.validTime(), eckit::mpi::comm());
  fset3d.shallowCopy(dx.increment().fieldSet());
  oops::FieldSet4D fset4d(fset3d);
  loc_->multiplySqrtAD(fset4d, dv.data(), 0);

  oops::Log::trace() << "Localization:multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << "Localization:multiply starting" << std::endl;

  // Create 4D fieldset
  oops::FieldSet4D fset4d({dx.validTime(), eckit::mpi::comm()});
  fset4d[0].shallowCopy(dx.increment().fieldSet());

  // SABER block chain multiplication
  loc_->multiply(fset4d);

  // ATLAS fieldset to Increment_
  dx.increment().synchronizeFields();

  oops::Log::trace() << "Localization:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  util::Timer timer(classname(), "print");
  os << "Localization:print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace saber
