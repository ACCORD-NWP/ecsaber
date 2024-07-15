/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <fstream>
#include <iostream>
#include <memory>

#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "eckit/utils/Translator.h"

#include "oops/runs/Test.h"

#include "util/Logger.h"

#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class GeometryTestFixture : TestFixture {
 public:
  GeometryTestFixture() {
    testconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "geometry"));
  }
  ~GeometryTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> testconf_;
};
// -----------------------------------------------------------------------------
CASE("test_geometry") {
  GeometryTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_geometry_constructor") {
    std::unique_ptr<quench::Geometry> resol(new quench::Geometry(*fix.testconf_));
    EXPECT(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_geometry_copy_constructor") {
    std::unique_ptr<quench::Geometry> xx(new quench::Geometry(*fix.testconf_));
    std::unique_ptr<quench::Geometry> resol(new quench::Geometry(*xx));
    EXPECT(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_geometry_methods") {
    std::unique_ptr<quench::Geometry> resol(new quench::Geometry(*fix.testconf_));
    EXPECT(resol->grid().uid() == fix.testconf_->getInt("uid"));
// TODO(Benjamin): add other tests
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
