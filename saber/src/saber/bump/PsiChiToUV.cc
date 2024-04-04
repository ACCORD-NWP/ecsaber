/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/PsiChiToUV.h"

#include "atlas/functionspace.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<PsiChiToUV> makerPsiChiToUV_("BUMP_PsiChiToUV");

// -----------------------------------------------------------------------------

namespace {

oops::patch::Variables createInnerVars(const oops::patch::Variables & outerVars) {
  oops::patch::Variables innerVars(std::vector<std::string>(
    {"stream_function", "velocity_potential"}));
  const int modelLevels(outerVars.getLevels("eastward_wind"));
  innerVars.addMetaData("stream_function", "levels", modelLevels);
  innerVars.addMetaData("velocity_potential", "levels", modelLevels);
  if (outerVars.hasMetaData("eastward_wind", "vert_coord")) {
    const std::string vert_coordName = outerVars.getMetaData<std::string>("eastward_wind",
      "vert_coord");
    innerVars.addMetaData("stream_function", "vert_coord", vert_coordName);
    innerVars.addMetaData("velocity_potential", "vert_coord", vert_coordName);
  }
  if (outerVars.hasMetaData("eastward_wind", "gmask")) {
    const std::string gmaskName = outerVars.getMetaData<std::string>("eastward_wind", "gmask");
    innerVars.addMetaData("stream_function", "gmask", gmaskName);
    innerVars.addMetaData("velocity_potential", "gmask", gmaskName);
  }
  return innerVars;
}

// -----------------------------------------------------------------------------

oops::patch::Variables createActiveVars(const oops::patch::Variables & innerVars,
                                 const oops::patch::Variables & outerVars) {
  oops::patch::Variables activeVars;
  activeVars += innerVars;
  activeVars += outerVars;
  const std::vector<std::string> activeStrings{"stream_function", "velocity_potential",
                                               "eastward_wind", "northward_wind"};
  activeVars.intersection(oops::patch::Variables(activeStrings));
  return activeVars;
}

}  // namespace

// -----------------------------------------------------------------------------

PsiChiToUV::PsiChiToUV(const oops::GeometryData & outerGeometryData,
                       const oops::patch::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(createInnerVars(outerVars)),
    outerVars_(outerVars),
    activeVars_(createActiveVars(innerVars_, outerVars_)),
    bumpParams_(params.calibrationParams.value() != boost::none ? *params.calibrationParams.value()
      : *params.readParams.value()),
    bump_(new BUMP(outerGeometryData, activeVars_, covarConf, bumpParams_,
      params.fieldsMetaData.value(), xb)) {
  oops::Log::trace() << classname() << "::PsiChiToUV starting" << std::endl;
  oops::Log::trace() << classname() << "::PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

PsiChiToUV::~PsiChiToUV() {
  oops::Log::trace() << classname() << "::~PsiChiToUV starting" << std::endl;
  util::Timer timer(classname(), "~PsiChiToUV");
  oops::Log::trace() << classname() << "::~PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyPsiChiToUV(fset);
  if (fset[innerVars_[0]].metadata().has("vert_coord")) {
    const std::string vert_coordName = fset[innerVars_[0]].metadata().getString("vert_coord");
    fset[outerVars_[0]].metadata().set("vert_coord", vert_coordName);
    fset[outerVars_[1]].metadata().set("vert_coord", vert_coordName);
  }
  if (fset[innerVars_[0]].metadata().has("gmask")) {
    const std::string gmaskName = fset[innerVars_[0]].metadata().getString("gmask");
    fset[outerVars_[0]].metadata().set("gmask", gmaskName);
    fset[outerVars_[1]].metadata().set("gmask", gmaskName);
  }
  fset.removeFields(innerVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(fset);
  if (fset[outerVars_[0]].metadata().has("vert_coord")) {
    const std::string vert_coordName = fset[outerVars_[0]].metadata().getString("vert_coord");
    fset[outerVars_[0]].metadata().set("vert_coord", vert_coordName);
    fset[outerVars_[1]].metadata().set("vert_coord", vert_coordName);
  }
  if (fset[outerVars_[0]].metadata().has("gmask")) {
    const std::string gmaskName = fset[outerVars_[0]].metadata().getString("gmask");
    fset[innerVars_[0]].metadata().set("gmask", gmaskName);
    fset[innerVars_[1]].metadata().set("gmask", gmaskName);
  }
  fset.removeFields(outerVars_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
