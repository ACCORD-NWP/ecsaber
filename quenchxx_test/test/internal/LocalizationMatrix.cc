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

#include "quench/LocalizationMatrix.h"
#include "quench/Geometry.h"
#include "quench/Variables.h"

#include "test/TestFixture.h"

namespace test {
class LocalizationMatrixFixture : TestFixture {
 public:
  LocalizationMatrixFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    vars_.reset(new quench::Variables(TestConfig::config()));
    eckit::LocalConfiguration cfg(TestConfig::config(), "Covariance");
    cfg_.reset(new eckit::LocalConfiguration(cfg));
  }
  ~LocalizationMatrixFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<quench::Variables> vars_;
  std::unique_ptr<const eckit::LocalConfiguration> cfg_;
};
// -----------------------------------------------------------------------------
CASE("test_LocalizationMatrix") {
  LocalizationMatrixFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_LocalizationMatrix_constructor") {
    std::unique_ptr<quench::LocalizationMatrix> locmat(
        new quench::LocalizationMatrix(*fix.resol_, *fix.vars_, *fix.cfg_));

    EXPECT(locmat.get() != NULL);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}

