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

#include "quench/Geometry.h"
#include "quench/GeoVaLs.h"
#include "quench/ObsOperator.h"
#include "quench/ObsSpace.h"
#include "quench/Variables.h"

#include "test/TestFixture.h"

#include "util/DateTime.h"

namespace test {

// -----------------------------------------------------------------------------
class GeoVaLsTestFixture : TestFixture {
 public:
  GeoVaLsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    ot_.reset(new quench::ObsSpace(otconf, *resol_, bgn, end));
    t1_.reset(new util::DateTime("2010-01-01T03:00:00Z"));
    t2_.reset(new util::DateTime("2010-01-02T06:00:00Z"));
    novar_.reset(new quench::Variables());
  }
  ~GeoVaLsTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<quench::ObsSpace> ot_;
  std::unique_ptr<util::DateTime> t1_;
  std::unique_ptr<util::DateTime> t2_;
  std::unique_ptr<quench::Variables> novar_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
CASE("test_GeoVaLs") {
  GeoVaLsTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_constructor") {
    std::unique_ptr<quench::GeoVaLs> gom(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    EXPECT(gom.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_nobs") {
    std::unique_ptr<quench::GeoVaLs> gom(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    EXPECT(gom->nobs() == 160);
  }
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_classname") {
    std::unique_ptr<quench::GeoVaLs> gom(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    EXPECT(gom->classname() == "quench::GeoVaLs");
  }
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_zero") {
    std::unique_ptr<quench::GeoVaLs> gom(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    gom->zero();
// TODO(Benjamin)
//    for (int i = 0; i < gom->nobs(); ++i) {
//      EXPECT((*gom)[i] == 0.0);
//    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_dot_product_with") {
    std::unique_ptr<quench::GeoVaLs> gom1(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    gom1->zero();
    std::unique_ptr<quench::GeoVaLs> gom2(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    gom2->zero();

    double zz = gom1->dot_product_with(*gom2);
    EXPECT(zz == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_GeoVaLs_fieldSet") {
    std::unique_ptr<quench::GeoVaLs> gom1(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    gom1->zero();
    std::unique_ptr<quench::GeoVaLs> gom2(new quench::GeoVaLs(*f.ot_, *f.novar_, *f.incr_, *f.t1_, *f.t2_));
    gom2->zero();
// TODO(Benjamin)
  }
// -----------------------------------------------------------------------------
}
// -----------------------------------------------------------------------------

}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
