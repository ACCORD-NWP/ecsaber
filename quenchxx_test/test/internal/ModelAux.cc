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
#include "eckit/utils/Translator.h"
#include "eckit/testing/Test.h"

#include "oops/runs/Test.h"

#include "util/Logger.h"

#include "quench/ModelAux.h"
#include "quench/ModelAuxCorrection.h"
#include "quench/ModelAuxCovariance.h"
#include "quench/Geometry.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasTestFixture : TestFixture {
 public:
  ModBiasTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelAux"));
    nobias_.reset(new eckit::LocalConfiguration());
    eckit::LocalConfiguration covarCfg(TestConfig::config(), "Covariance");
    quench::ModelAuxCovariance covar(covarCfg, *resol_);
    bias1_ = biasconf_->getDouble("bias");
    bias2_ = 2.5 * bias1_;
    dbias_.reset(new quench::ModelAuxCorrection(*resol_, covar.config()));
    dbias_->bias() = bias2_;
  }
  ~ModBiasTestFixture() {}
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> biasconf_;
  std::unique_ptr<const eckit::LocalConfiguration> nobias_;
  std::unique_ptr<quench::ModelAuxCorrection> dbias_;
  double bias1_;
  double bias2_;
};
// -----------------------------------------------------------------------------
CASE("test_ModelAux") {
  ModBiasTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_bias") {
    quench::ModelAux mbias(*fix.resol_, *fix.biasconf_);
    EXPECT(mbias.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_nobias") {
    quench::ModelAux mbias(*fix.resol_, *fix.nobias_);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias.bias() == 0.0);
}
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_active") {
    quench::ModelAux mbias1(*fix.resol_, *fix.biasconf_);
    quench::ModelAux mbias2(*fix.resol_, mbias1);
    EXPECT(mbias2.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_inactive") {
    quench::ModelAux mbias1(*fix.resol_, *fix.nobias_);
    quench::ModelAux mbias2(*fix.resol_, mbias1);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_active_copy") {
    quench::ModelAux mbias1(*fix.resol_, *fix.biasconf_);
    quench::ModelAux mbias2(mbias1, true);
    EXPECT(mbias2.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_active_nocopy") {
    quench::ModelAux mbias1(*fix.resol_, *fix.biasconf_);
    quench::ModelAux mbias2(mbias1, false);

    // because copy flag is set to false, bias is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_inactive_copy") {
    quench::ModelAux mbias1(*fix.resol_, *fix.nobias_);
    quench::ModelAux mbias2(mbias1, true);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_constructor_inactive_nocopy") {
    quench::ModelAux mbias1(*fix.resol_, *fix.nobias_);
    quench::ModelAux mbias2(mbias1, false);
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_classname") {
    quench::ModelAux mbias(*fix.resol_, *fix.biasconf_);
    EXPECT(mbias.classname() == "quench::ModelAux");
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_compound_assignment_active") {
    quench::ModelAux mbias(*fix.resol_, *fix.biasconf_);

    mbias += *fix.dbias_;

    EXPECT(mbias.bias() == fix.bias1_ + fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_compound_assignment_inactive") {
    quench::ModelAux mbias(*fix.resol_, *fix.nobias_);

    mbias += *fix.dbias_;

    // because the biasconf_ is empty,
    // mbias bias value is not modified from 0.0
    EXPECT(mbias.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_bias") {
    quench::ModelAux mbias(*fix.resol_, *fix.biasconf_);
    EXPECT(mbias.bias() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ModelAux_stream_output") {
    quench::ModelAux mbias(*fix.resol_, *fix.biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelAuxTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << mbias;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = fix.bias1_;
    double bias = 0.0;
    int biasStartPos = 12;  // length of "ModelAux = " is 12
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

      EXPECT(oops::is_close(testBias, bias, 0.0001));
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
