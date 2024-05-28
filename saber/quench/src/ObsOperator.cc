/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ObsOperator.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "util/DateTime.h"
#include "util/Logger.h"

#include "src/GeoVaLs.h"
#include "src/ObsAuxControl.h"
#include "src/ObsVec.h"
#include "src/Traits.h"
#include "src/Variables.h"

namespace quench {

// -----------------------------------------------------------------------------

static oops::ObsOperatorMaker<Traits, ObsOperator> makerObsOperatorDefault_("default");

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(const ObsSpace & obsSpace,
                         const eckit::Configuration & conf)
  : obsSpace_(obsSpace), inputs_(new Variables(conf)) {
  oops::Log::trace() << classname() << "::obsEquiv" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsOperator::obsEquiv(const GeoVaLs & gv, ObsVec & ovec,
                           const ObsAuxControlPtrMap_ & bias) const {
  oops::Log::trace() << classname() << "::obsEquiv starting" << std::endl;

  // Check number of GeoVaLs variables
  ASSERT(gv.fieldSet().size() == 1);

  // Get GeoVaLs view
  const auto gvField = gv.fieldSet()[0];
  const auto gvView = atlas::array::make_view<double, 1>(gvField);

  // Get bias
  double bias_ = 0.0;  // TODO(Benjamin): bias should also be an atlas fieldset
  using icst_ = typename ObsAuxControlPtrMap_::const_iterator;
  icst_ it = bias.find("ObsAuxControl");
  if (it != bias.end()) {
    const ObsAuxControl * pbias = dynamic_cast<const ObsAuxControl*>(it->second.get());
    ASSERT(pbias != nullptr);
    bias_ = (*pbias).value();
  }

  // Compute observation equivalent
  for (int jo = 0; jo < gvField.shape(0); ++jo) {
    const int ii = gv.obsIndex(jo);
    ovec(ii) = gvView(jo)+bias_;
  }

  oops::Log::trace() << classname() << "::obsEquiv done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "ObsOperator: quench";

  oops::Log::trace() << classname() << "::print end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
