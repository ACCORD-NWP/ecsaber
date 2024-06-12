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
#include <memory> // for std::unique_ptr

#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "eckit/utils/Translator.h"

#include "oops/runs/Test.h"

#include "util/Logger.h"

#include "quench/ObsAux.h"
#include "quench/ObsAuxIncrement.h"
#include "quench/ObsSpace.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsAuxTestFixture : TestFixture {
 public:
  ObsAuxTestFixture() {
    off_.reset(new eckit::LocalConfiguration());
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsauxCovariance"));
    eckit::LocalConfiguration bconf(TestConfig::config(), "ObsAux");
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    obstable_.reset(new quench::ObsSpace(otconf, *resol_, bgn, end));
    obias_.reset(new quench::ObsAux(*obstable_, bconf));
    bias1_ = bconf.getDouble("bias");
    bias2_ = 3.5 * bias1_;
    fact_ = 1.234;
  }
  ~ObsAuxTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> off_;
  std::unique_ptr<const eckit::LocalConfiguration> conf_;
  double bias1_;
  double bias2_;
  double fact_;
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<quench::ObsAux> obias_;
  std::unique_ptr<quench::ObsSpace> obstable_;
};

// -----------------------------------------------------------------------------
CASE("test_ObsAuxIncrement") {
  ObsAuxTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_constructor_active") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.conf_);
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_constructor_inactive") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.off_);
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_ctor_active_copy") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, true);

    EXPECT(dob2.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_ctor_active_no_copy") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, false);

    // because the copy is false,
    // the active_ flag is true and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_ctor_inactive_copy") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_ctor_inactive_no_copy") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_constructor_config") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, *fix.conf_);

    EXPECT(dob2.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_copy_constructor_no_config") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;

    quench::ObsAuxIncrement dob2(dob1, eckit::LocalConfiguration());

    // because the covarCfg is empty when used (regardless of the cfg),
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_diff_active") {
    quench::ObsAux obias2(*fix.obstable_, *fix.conf_);
    obias2.value() = fix.bias2_;

    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.conf_);

    dob.diff(*fix.obias_, obias2);

    EXPECT(dob.value() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_diff_inactive") {
    quench::ObsAux obias2(*fix.obstable_, *fix.conf_);
    obias2.value() = fix.bias2_;

    // construct the dob object with empty config
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.off_);

    dob.diff(*fix.obias_, obias2);

    // because the OBC has empty config the active_flag is false and
    // the diff will not be performed, leaving the OBC value unchanged
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_zero") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.conf_);
    dob.value() = fix.bias1_;

    dob.zero();

    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_active") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 = dob2;

    EXPECT(dob1.value() == fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_inactive") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 = dob2;

    // because the OBC has empty config the active_ flag is false
    EXPECT(dob1.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_add_active") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 += dob2;

    EXPECT(dob1.value() == fix.bias1_ + fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assign_add_inactive") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 += dob2;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_subtract_active") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 -= dob2;

    EXPECT(dob1.value() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_subtract_inactive") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 -= dob2;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_multiply_active") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.conf_);
    dob.value() = fix.bias1_;

    dob *= fix.fact_;

    EXPECT(dob.value() == fix.bias1_ * fix.fact_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_assignment_multiply_inactive") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.off_);
    dob.value() = fix.bias1_;

    dob *= fix.fact_;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_axpy_active") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1.axpy(fix.fact_, dob2);

    EXPECT(dob1.value() == fix.bias1_ + (fix.fact_ * fix.bias2_));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_axpy_inactive") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1.axpy(fix.fact_, dob2);

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_dot_product_with_active") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.conf_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    EXPECT(dpwResult == fix.bias1_ * fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_dot_product_with_inactive") {
    quench::ObsAuxIncrement dob1(*fix.obstable_, *fix.off_);
    dob1.value() = fix.bias1_;
    quench::ObsAuxIncrement dob2(dob1);
    dob2.value() = fix.bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    // because the OBC has empty config, the result is 0
    EXPECT(dpwResult == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAuxIncrement_stream_output") {
    quench::ObsAuxIncrement dob(*fix.obstable_, *fix.conf_);
    dob.value() = fix.bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsAuxIncrementTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double bias = 0.0;
    int biasStartPos = 20;  // length of "ObsAuxIncrement = " is 20
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);  // ignore first (blank) line
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = eckit::Translator<std::string, double>()(inputBias);
      }
      catch(eckit::BadParameter const&) {
        oops::Log::error() << "operator<< incorrectly output a non-double" << std::endl;
      }

      EXPECT(oops::is_close(fix.bias1_, bias, 0.0001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
