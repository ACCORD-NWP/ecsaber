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

#include "quench/ModelAux.h"
#include "quench/ModelAuxCorrection.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasTestFixture : TestFixture {
 public:
  ModBiasTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelAuxCovariance"));
    nobias_.reset(new eckit::LocalConfiguration());
    bias1_ = TestConfig::config().getDouble("ModelAux.bias");
    bias2_ = 2.5 * bias1_;
    fact_ = 1.2345;
  }
  ~ModBiasTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> conf_;
  std::unique_ptr<const eckit::LocalConfiguration> nobias_;
  double bias1_;
  double bias2_;
  double fact_;
};
// -----------------------------------------------------------------------------
CASE("test_ModelAuxCorrection") {
  ModBiasTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_constructor_config") {
    std::unique_ptr<quench::ModelAuxCorrection> dx(
      new quench::ModelAuxCorrection(*fix.resol_, *fix.conf_));
    EXPECT(dx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_constructor_no_config") {
    std::unique_ptr<quench::ModelAuxCorrection> dx(
      new quench::ModelAuxCorrection(*fix.resol_, *fix.nobias_));
    EXPECT(dx.get() != NULL);
}
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_active_copy") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;

    quench::ModelAuxCorrection dx2(dx1, true);

    EXPECT(dx2.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_active_no_copy") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;

    // construct a copy of it with the copy flag set to false
    quench::ModelAuxCorrection dx2(dx1, false);

    // because the copy is false,
    // the active_ flag is true and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_inactive_copy") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;

    // construct a copy of it with the copy flag set to true
    quench::ModelAuxCorrection dx2(dx1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_inactive_no_copy") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;

    // construct a copy of it with the copy flag set to false
    quench::ModelAuxCorrection dx2(dx1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_config_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;

    quench::ModelAuxCorrection dx2(dx1, *fix.conf_);

    EXPECT(dx2.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_copy_ctor_config_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;

    quench::ModelAuxCorrection dx2(dx1, *fix.conf_);

     // because the covarCfg is empty when used (regardless of the cfg),
     // the active_ flag is false and the bias_ value is 0.0
     EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_diff_active") {
    quench::ModelAux xx1(*fix.resol_, *fix.conf_);
    xx1.bias() = fix.bias1_;
    quench::ModelAux xx2(*fix.resol_, *fix.conf_);
    xx2.bias() = fix.bias2_;

    quench::ModelAuxCorrection dx(*fix.resol_, *fix.conf_);

    dx.diff(xx1, xx2);

    EXPECT(dx.bias() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_diff_inactive") {
    quench::ModelAux xx1(*fix.resol_, *fix.conf_);
    xx1.bias() = fix.bias1_;
    quench::ModelAux xx2(*fix.resol_, *fix.conf_);
    xx2.bias() = fix.bias2_;

    quench::ModelAuxCorrection dx(*fix.resol_, *fix.nobias_);

    dx.diff(xx1, xx2);

    // because the active_ flag is false, the bias cannot be updated
    EXPECT(dx.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_zero") {
    quench::ModelAuxCorrection dx(*fix.resol_, *fix.conf_);
    dx.bias() = fix.bias1_;

    dx.zero();

    EXPECT(dx.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    dx1 = dx2;

    // the original MBC should have the same bias value as the copy MBC
    EXPECT(dx1.bias() == fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;

    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    dx1 = dx2;

    // the active_ value is zero, so the bias will be zero
    EXPECT(dx1.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_add_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    dx1 += dx2;

    EXPECT(dx1.bias() == fix.bias1_ + fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_add_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.nobias_);
    dx2.bias() = fix.bias2_;

    dx1 += dx2;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_subtract_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    dx1 -= dx2;

    EXPECT(dx1.bias() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_subtract_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.nobias_);
    dx2.bias() = fix.bias2_;

    dx1 -= dx2;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_multiply_active") {
    quench::ModelAuxCorrection dx(*fix.resol_, *fix.conf_);
    dx.bias() = fix.bias1_;

    dx *= fix.fact_;

    EXPECT(dx.bias() == fix.bias1_ * fix.fact_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_assignment_multiply_inactive") {
    quench::ModelAuxCorrection dx(*fix.resol_, *fix.nobias_);
    dx.bias() = fix.bias1_;

    dx *= fix.fact_;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_axpy_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    dx1.axpy(fix.fact_, dx2);

    EXPECT(dx1.bias() == (fix.bias1_ + fix.fact_ * fix.bias2_));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_axpy_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.nobias_);
    dx2.bias() = fix.bias2_;

    dx1.axpy(fix.fact_, dx2);

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_dot_product_with_active") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.conf_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.conf_);
    dx2.bias() = fix.bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    EXPECT(dpwResult == fix.bias1_ * fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_dot_product_with_inactive") {
    quench::ModelAuxCorrection dx1(*fix.resol_, *fix.nobias_);
    dx1.bias() = fix.bias1_;
    quench::ModelAuxCorrection dx2(*fix.resol_, *fix.nobias_);
    dx2.bias() = fix.bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    // because of the empty config, the result is 0.0
    EXPECT(dpwResult == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_read") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_write") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_stream_output") {
    quench::ModelAuxCorrection dx(*fix.resol_, *fix.conf_);
    dx.bias() = fix.bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelAuxCorrectionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dx;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = fix.bias1_;
    double bias = 0.0;
    int biasStartPos = 22;  // length of "ModelAuxCorrection = " is 22
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

      EXPECT(oops::is_close(testBias, bias, 0.000001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAuxCorrection_bias") {
    quench::ModelAuxCorrection dx(*fix.resol_, *fix.conf_);
    dx.bias() = fix.bias1_;

    // this one test checks both the setting and getting of the bias value
    EXPECT(dx.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
