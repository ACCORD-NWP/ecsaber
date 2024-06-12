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

#include "quench/ObsAuxIncrement.h"
#include "quench/ObsauxCovariance.h"
#include "quench/ObsSpace.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsAuxTestFixture : TestFixture {
 public:
  ObsAuxTestFixture() {
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsAux"));
    nobias_.reset(new eckit::LocalConfiguration());
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    obstable_.reset(new quench::ObsSpace(otconf, *resol_, bgn, end));
  }
  ~ObsAuxTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> biasconf_;
  std::unique_ptr<const eckit::LocalConfiguration> nobias_;
  std::unique_ptr<const eckit::LocalConfiguration> covconf_;
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<quench::ObsSpace> obstable_;
};
// -----------------------------------------------------------------------------
CASE("test_ObsauxCovariance") {
  ObsAuxTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_constructor_conf") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.covconf_);
    EXPECT(obcovar.active() == true);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_constructor_no_conf") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.nobias_);
    EXPECT(obcovar.active() == false);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_linearize") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_multiply_active") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.covconf_);

    quench::ObsAuxIncrement db1(*fix.obstable_, *fix.covconf_);
    db1.value() = 2.0;
    quench::ObsAuxIncrement db2(db1, *fix.covconf_);

    obcovar.multiply(db1, db2);

    const double stdev = fix.covconf_->getDouble("standard_deviation");
    EXPECT(db2.value() == db1.value() * stdev * stdev);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_multiply_inactive") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.nobias_);

    quench::ObsAuxIncrement db1(*fix.obstable_, *fix.nobias_);
    db1.value() = 2.0;
    quench::ObsAuxIncrement db2(db1, *fix.covconf_);

    obcovar.multiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    EXPECT(db2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_invMult_active") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.covconf_);

    quench::ObsAuxIncrement db1(*fix.obstable_, *fix.covconf_);
    db1.value() = 2.0;
    quench::ObsAuxIncrement db2(db1, *fix.covconf_);

    obcovar.inverseMultiply(db1, db2);

    const double stdev = fix.covconf_->getDouble("standard_deviation");
    EXPECT(db2.value() == db1.value() * 1.0 / (stdev * stdev));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsauxCovariance_invMult_inactive") {
    quench::ObsauxCovariance obcovar(*fix.obstable_, *fix.nobias_);

    quench::ObsAuxIncrement db1(*fix.obstable_, *fix.nobias_);
    db1.value() = 2.0;
    quench::ObsAuxIncrement db2(db1, *fix.covconf_);

    obcovar.inverseMultiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    EXPECT(db2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
