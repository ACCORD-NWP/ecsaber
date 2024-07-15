/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>

#include <memory> // for std::unique_ptr

#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/runs/Test.h"

#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

#include "quench/FieldL95.h"
#include "quench/ModelAux.h"
#include "quench/Model.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModelTestFixture : TestFixture {
 public:
  ModelTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    nlconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "model"));
  }
  ~ModelTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> nlconf_;
};
// -----------------------------------------------------------------------------
CASE("test_Model") {
ModelTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_Model_constructor") {
    std::unique_ptr<quench::Model> model(new quench::Model(*fix.resol_, *fix.nlconf_));
    EXPECT(model != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Model_get_classname") {
    quench::Model model(*fix.resol_, *fix.nlconf_);
    EXPECT(model.classname() == "quench::Model");
  }
// -----------------------------------------------------------------------------
  SECTION("test_Model_get_timestep") {
    quench::Model model(*fix.resol_, *fix.nlconf_);
    util::Duration dt(fix.nlconf_->getString("tstep"));
    EXPECT(model.timeGeometry().toSeconds() == dt.toSeconds());
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
