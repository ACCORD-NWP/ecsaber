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

#include "quench/ModelAuxCorrection.h"
#include "quench/ModelAuxCovariance.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasCovTestFixture : TestFixture {
 public:
  ModBiasCovTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelAuxCovariance"));
    nobias_.reset(new eckit::LocalConfiguration());
  }
  ~ModBiasCovTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> covconf_;
  std::unique_ptr<const eckit::LocalConfiguration> nobias_;
};
// -----------------------------------------------------------------------------
CASE("test_ModelAuxCovariance") {
  ModBiasCovTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_constructor_conf") {
    quench::ModelAuxCovariance bcovar(*f.covconf_, *f.resol_);
    EXPECT(bcovar.active() == true);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_constructor_no_conf") {
    quench::ModelAuxCovariance bcovar(*f.nobias_, *f.resol_);
    EXPECT(bcovar.active() == false);
}
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_linearize") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_multiply_activei") {
    // construct the ModelAuxCorrection object
    quench::ModelAuxCovariance bcovar(*f.covconf_, *f.resol_);
    quench::ModelAuxCorrection dbias1(*f.resol_, *f.covconf_);
    dbias1.bias() = 2.0;
    quench::ModelAuxCorrection dbias2(dbias1, true);

    bcovar.multiply(dbias1, dbias2);

    double stdev = f.covconf_->getDouble("standard_deviation");
    EXPECT(dbias2.bias() == dbias1.bias() * stdev * stdev);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_multiply_inactive") {
    // construct the ModelAuxCorrection object
    quench::ModelAuxCovariance bcovar(*f.nobias_, *f.resol_);
    quench::ModelAuxCorrection dbias1(*f.resol_, *f.covconf_);
    dbias1.bias() = 2.0;
    quench::ModelAuxCorrection dbias2(dbias1, true);

    bcovar.multiply(dbias1, dbias2);

    // because the covconf_ is empty, the bias is set to 0
    EXPECT(dbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_invMult_active") {
    // construct the ModelAuxCorrection object
    quench::ModelAuxCovariance bcovar(*f.covconf_, *f.resol_);
    quench::ModelAuxCorrection dbias1(*f.resol_, *f.covconf_);
    dbias1.bias() = 2.0;
    quench::ModelAuxCorrection dbias2(dbias1, true);

    bcovar.inverseMultiply(dbias1, dbias2);

    double stdev = f.covconf_->getDouble("standard_deviation");
    EXPECT(dbias2.bias() == dbias1.bias() *1.0 / (stdev * stdev));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_invMult_inactive") {
    // construct the ModelAuxCorrection object
    quench::ModelAuxCovariance bcovar(*f.nobias_, *f.resol_);
    quench::ModelAuxCorrection dbias1(*f.resol_, *f.covconf_);
    dbias1.bias() = 2.0;
    quench::ModelAuxCorrection dbias2(dbias1, true);

    bcovar.inverseMultiply(dbias1, dbias2);

    // because the covconf_ is empty, the bias is set to 0
    EXPECT(dbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCovariance_active") {
    quench::ModelAuxCovariance bcovar(*f.covconf_, *f.resol_);
    EXPECT(bcovar.active() == true);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}

