/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory> // for std::unique_ptr

#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/runs/Test.h"

#include "quench/ObsOperator.h"
#include "quench/ObsSpace.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsTestFixture : TestFixture {
 public:
  ObsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    ot_.reset(new quench::ObsSpace(otconf, *resol_, bgn, end));
  }
  ~ObsTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<quench::ObsSpace> ot_;
};
// -----------------------------------------------------------------------------
CASE("test_ObsOperator") {
  ObsTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsOperator_constructor") {
    std::unique_ptr<quench::ObsOperator>
      obs(new quench::ObsOperator(*fix.ot_, TestConfig::config()));
    EXPECT(obs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_classname") {
    std::unique_ptr<quench::ObsOperator>
      obs(new quench::ObsOperator(*fix.ot_, TestConfig::config()));
    EXPECT(obs->classname() == "quench::ObsOperator");
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
