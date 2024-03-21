/*
 * (C) Copyright 2024 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>

#include "saber/generic/FakeLevels.h"

#include "eckit/config/LocalConfiguration.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<FakeLevels> makerFakeLevels_("fake levels");

// -----------------------------------------------------------------------------

FakeLevels::FakeLevels(const oops::GeometryData & outerGeometryData,
                       const oops::patch::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    params_(params),
    nz_(params.fakeLevels.value().size()),
    validTime_(xb.validTime()),
    gdata_(outerGeometryData),
    comm_(gdata_.comm()),
    outerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars_)),
    suffix_("_fakeLevels"),
    innerVars_(createInnerVars(nz_, activeVars_, outerVars)),
    fakeLevels_(params_.fakeLevels.value()),
    fieldsMetaData_(params_.fieldsMetaData.value()) {
  oops::Log::trace() << classname() << "::FakeLevels starting" << std::endl;
  oops::Log::trace() << classname() << "::FakeLevels done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet outerFset;

  for (const auto & outerVar : activeVars_.variables()) {
    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Get inner field
    const auto innerView = atlas::array::make_view<double, 2>(fset[innerVar]);

    // Create outer field
    atlas::Field outerField = gdata_.functionSpace().createField<double>(
      atlas::option::name(outerVar) | atlas::option::levels(1) | atlas::option::halo(1));
    auto outerView = atlas::array::make_view<double, 2>(outerField);

    // Get weight
    const auto weightsView = atlas::array::make_view<double, 2>((*weights_)[outerVar]);

    // Reduce to a single level
    outerView.assign(0.0);
    for (int jnode = 0; jnode < outerField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          outerView(jnode, 0) += innerView(jnode, k)*weightsView(jnode, k);
        }
      }
    }

    // Add field
    outerFset.add(outerField);
  }

  // Copy outer fieldset
  fset.fieldSet() = outerFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet innerFset;

  for (const auto & outerVar : activeVars_.variables()) {
    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Get outer field
    const auto outerView = atlas::array::make_view<double, 2>(fset[outerVar]);
    ASSERT(fset[outerVar].shape(1) == 1);

    // Create inner field
    atlas::Field innerField = gdata_.functionSpace().createField<double>(
      atlas::option::name(innerVar) | atlas::option::levels(nz_) | atlas::option::halo(1));
    auto innerView = atlas::array::make_view<double, 2>(innerField);

    // Get weight
    const auto weightsView = atlas::array::make_view<double, 2>((*weights_)[outerVar]);

    // Extend to multiple levels
    innerView.assign(0.0);
    for (int jnode = 0; jnode < innerField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          innerView(jnode, k) = outerView(jnode, 0)*weightsView(jnode, k);
        }
      }
    }

    // Add field
    innerFset.add(innerField);
  }

  // Copy inner fieldset
  fset.fieldSet() = innerFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> FakeLevels::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (params_.inputModelFilesConf.value() != boost::none) {
    for (const auto & conf : *params_.inputModelFilesConf.value()) {
      // File name
      const std::string fileName = conf.getString("parameter");

      // Add pair
      inputs.push_back(std::make_pair(fileName, conf));
    }
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void FakeLevels::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create rv
  oops::FieldSet3D rv(validTime_, comm_);

  // Get rv from input files if present
  for (size_t ji = 0; ji < fsetVec.size(); ++ji) {
    if (fsetVec[ji].name() == "rv") {
      for (const auto & var : activeVars_.variables()) {
        // Save field
        rv.add(fsetVec[ji][var]);
      }
    }
  }

  // Get vertical support from yaml if needed
  if (rv.empty()) {
    ASSERT(params_.rvFromYaml.value() != boost::none);
    for (const auto & var : activeVars_.variables()) {
      atlas::Field rvField = gdata_.functionSpace().createField<double>(
        atlas::option::name(var) | atlas::option::levels(1));
      auto rvView = atlas::array::make_view<double, 2>(rvField);
      for (int jnode = 0; jnode < rvField.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          rvView(jnode, 0) = *params_.rvFromYaml.value();
        }
      }
      rv.add(rvField);
    }
  }

  // Prepare weights
  weights_.reset(new oops::FieldSet3D(validTime_, comm_));
  for (const auto & var : activeVars_.variables()) {
    // Check number of levels of outer variables
    ASSERT(outerVars_.getLevels(var) == 1);

    // Get vertical coordinate
    const std::string key = var + ".vert_coord";
    const std::string vertCoordName = fieldsMetaData_.getString(key, "vert_coord");
    const atlas::Field vertCoordField = gdata_.fieldSet()[vertCoordName];
    const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoordField);

    // Get vertical support
    const auto rvView = atlas::array::make_view<double, 2>(rv[var]);

    // Create weight field
    atlas::Field field = gdata_.functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(nz_) | atlas::option::halo(1));

    // Set to zero
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);

    // Compute weights
    std::vector<double> wgt(nz_);
    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Check fake levels extrema
        ASSERT(vertCoordView(jnode, 0) >= fakeLevels_[0]);
        ASSERT(vertCoordView(jnode, 0) <= fakeLevels_[nz_-1]);

        // Compute raw weights
        double wgtSum = 0.0;
        for (size_t k = 0; k < nz_; ++k) {
          const double normDist = std::abs(vertCoordView(jnode, 0)-fakeLevels_[k])/rvView(jnode, 0);
          if (fakeLevels_[k] < vertCoordView(jnode, 0)) {
            // Under ground
            wgt[k] = 0.0;
          } else {
            // Above ground
            wgt[k] = oops::gc99(normDist);
            wgt[k] = std::max(0.0, wgt[k]);
          }
          wgtSum += wgt[k];
        }

        // Normalize weights
        if (wgtSum > 0) {
          for (size_t k = 0; k < nz_; ++k) {
            wgt[k] /= wgtSum;
          }
        } else {
          throw eckit::UserError("fake levels are too far apart", Here());
        }

        // Weight square-roots
        for (size_t k = 0; k < nz_; ++k) {
          ASSERT(wgt[k] >= 0.0);
          wgt[k] = std::sqrt(wgt[k]);
        }

        // Copy weights
        for (size_t k = 0; k < nz_; ++k) {
          view(jnode, k) = wgt[k];
        }
      }
    }

    // Halo exchange
    field.haloExchange();

    // Add field
    weights_->add(field);
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::patch::Variables FakeLevels::createInnerVars(const atlas::idx_t & innerVerticalLevels,
                                            const oops::patch::Variables & activeVars,
                                            const oops::patch::Variables & outerVars) const {
  oops::Log::trace() << classname() << "::createInnerVars starting" << std::endl;

  oops::patch::Variables innerVars;
  for (const std::string & varName : outerVars.variables()) {
    if (activeVars.has(varName)) {
      const std::string newVarName = varName + suffix_;
      innerVars.push_back(newVarName);
      innerVars.addMetaData(newVarName, "levels", innerVerticalLevels);
    } else {
      innerVars.push_back(varName);
    }
  }

  oops::Log::trace() << classname() << "::createInnerVars done" << std::endl;
  return innerVars;
}

// -----------------------------------------------------------------------------

void FakeLevels::print(std::ostream & os) const {
  os << classname();
}

}  // namespace generic
}  // namespace saber
