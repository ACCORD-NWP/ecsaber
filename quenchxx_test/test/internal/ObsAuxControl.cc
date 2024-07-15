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
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsAux"));
    nobias_.reset(new eckit::LocalConfiguration());
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
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
CASE("test_ObsAux") {
  ObsAuxTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_constructor_bias") {
    quench::ObsAux ob(*fix.obstable_, *fix.biasconf_);
    EXPECT(ob.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_constructor_no_bias") {
    quench::ObsAux ob(*fix.obstable_, *fix.nobias_);
    EXPECT(ob.value() == 0.0);
}
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_copy_constructor_config_copy") {
    quench::ObsAux ob1(*fix.obstable_, *fix.biasconf_);
    quench::ObsAux ob2(ob1, true);
    EXPECT(ob2.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_copy_constructor_config_no_copy") {
    quench::ObsAux ob1(*fix.obstable_, *fix.biasconf_);
    quench::ObsAux ob2(ob1, false);

    // bias value is 0 because copy flag was set to false
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_copy_constructor_no_config_copy") {
    quench::ObsAux ob1(*fix.obstable_, *fix.nobias_);
    quench::ObsAux ob2(ob1, true);

    // bias value is 0 because an empty config was used
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_copy_constructor_no_config_no_copy") {
    quench::ObsAux ob1(*fix.obstable_, *fix.nobias_);
    quench::ObsAux ob2(ob1, true);

    // bias value is 0 because an empty config was used
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_compound_assignment_add_active") {
    quench::ObsAux ob(*fix.obstable_, *fix.biasconf_);

    // construct an ObsAuxIncrement object
    quench::ObsAuxIncrement ObsAuxIncrement(*fix.obstable_, *fix.covconf_);
    ObsAuxIncrement.value() = 3.14;

    ob += ObsAuxIncrement;

    EXPECT(ob.value() == fix.biasconf_->getDouble("bias") + 3.14);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_compound_assignment_add_inactive") {
    quench::ObsAux ob(*fix.obstable_, *fix.nobias_);

    // construct an ObsAuxIncrement object
    quench::ObsAuxIncrement ObsAuxIncrement(*fix.obstable_, *fix.covconf_);
    ObsAuxIncrement.value() = 3.14;

    ob += ObsAuxIncrement;

    // the bias value will be 0 because the ob had an empty config
    EXPECT(ob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_value") {
    quench::ObsAux ob(*fix.obstable_, *fix.biasconf_);
    EXPECT(ob.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsAux_stream_output") {
    quench::ObsAux ob(*fix.obstable_, *fix.biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsAuxTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << ob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = fix.biasconf_->getDouble("bias");
    double bias = 0.0;
    int biasStartPos = 10;  // length of "ObsAux = " is 10
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);
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
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
