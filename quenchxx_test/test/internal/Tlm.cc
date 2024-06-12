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

#include "oops/runs/Test.h"

#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

#include "quench/ModelAux.h"
#include "quench/ModelAuxCorrection.h"
#include "quench/ModelAuxCovariance.h"
#include "quench/Model.h"
#include "quench/ModelTrajectory.h"
#include "quench/Geometry.h"
#include "quench/State.h"
#include "quench/Tlm.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class TlmTestFixture : TestFixture {
 public:
  TlmTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    eckit::LocalConfiguration mod(TestConfig::config(), "model");
    model_.reset(new quench::Model(*resol_, mod));
    tlconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "linearmodel"));
  }
  ~TlmTestFixture() {}
  std::unique_ptr<quench::Model>   model_;
  std::unique_ptr<quench::Geometry> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> tlconf_;
};
// -----------------------------------------------------------------------------
CASE("test_Tlm") {
  TlmTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_constructor") {
    std::unique_ptr<quench::Tlm> tlm(new quench::Tlm(*fix.resol_, *fix.tlconf_));
    EXPECT(tlm.get() != NULL);
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_Tlm_set_get_Trajectory") {
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_, *fix.model_);

    // construct a State object
    std::string date_string("2014-09-12T09:35:00Z");
    util::DateTime dt(date_string);
    oops::Variables vars();
    std::unique_ptr<quench::State> State(new quench::State(*fix.resol_, vars, dt));

    // construct a ModelAux object
    std::shared_ptr<const eckit::LocalConfiguration> biasCfg(new eckit::LocalConfiguration(TestConfig::config(), "ModelAux"));
    std::unique_ptr<quench::ModelAux> ModelAux(new quench::ModelAux(*biasCfg, *fix.resol_));

    tlm.setTrajectory(*State, *ModelAux);

    std::unique_ptr<quench::ModelTrajectory> modelTraj(new quench::ModelTrajectory(true));
    oops::Log::error() << "PMC: getTraj result <" << tlm.getTrajectory(dt)->get(1) << ">" << std::endl;
    modelTraj->set(tlm.getTrajectory(dt)->get(1));
    //oops::Log::error() << "PMC: modelTraj 1st " << modelTraj->get(1)
    oops::Log::error() << "PMC: modelTraj resol  <" << modelTraj->get(1).resol() << ">" << std::endl;
    for(int i = 0; i < modelTraj->get(1).resol(); ++i) {
      oops::Log::error() << "PMC: modelTraj FieldL95 elements " << modelTraj->get(1)[i] << std::endl;
    }
  }
*/
// -----------------------------------------------------------------------------
/*
  SECTION("test_Tlm_stepTL") {
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_, *fix.model_);

    // construct a FieldL95 object
    std::unique_ptr<quench::FieldL95> fieldL95(new quench::FieldL95(*fix.resol_));

    // construct a datetime object
    std::string dateString("2014-10-08T11:25:45Z");
    util::DateTime dt(dateString);

    //construct a ModelAuxCorrection object
    std::unique_ptr<const eckit::LocalConfiguration> covarCfg(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
    std::unique_ptr<quench::ModelAuxCovariance> ModelAuxCovariance(
        new quench::ModelAuxCovariance(*covarCfg, *fix.resol_));
    std::unique_ptr<quench::ModelAuxCorrection> ModelAuxCorrection(
        new quench::ModelAuxCorrection(*ModelAuxCovariance));

    tlm.stepTL(*fieldL95, dt, *ModelAuxCorrection);

    oops::Log::error() << "PMC: original dt was " << dateString << std::endl;
    oops::Log::error() << "PMC: modified dt is  " << dt.toString() << std::endl;
  }
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_stepAD") {
  }
*/
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_get_classname") {
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_);
    EXPECT(tlm.classname() == "quench::Tlm");
  }
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_get_timestep") {
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_);
    util::Duration dt(fix.tlconf_->getString("tstep"));
    EXPECT(tlm.timeGeometry().toSeconds() == dt.toSeconds());
  }
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_get_geometry") {
    eckit::LocalConfiguration rescf(TestConfig::config(), "geometry");
    quench::Geometry resol(rescf);
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_);
    EXPECT(tlm.geometry().npoints() == resol.npoints());
  }
// -----------------------------------------------------------------------------
  SECTION("test_Tlm_stream_output") {
    quench::Tlm tlm(*fix.resol_, *fix.tlconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("TlmTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << tlm;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputStartsWith;
    std::string testString("Tlm: resol = ");
    int endPos = testString.size();
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);

      inputStartsWith = inputString.substr(0, endPos);

      EXPECT(inputStartsWith == testString);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() <<"operator<< functionality cannot be determined" << std::endl;
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
