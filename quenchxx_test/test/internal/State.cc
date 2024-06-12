/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory> // for std::unique_ptr
#include <sstream>

#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/runs/Test.h"

#include "util/DateTime.h"
#include "util/Logger.h"

#include "quench/FieldL95.h"
#include "quench/GeoVaLs.h"
#include "quench/Increment.h"
#include "quench/Locations.h"
#include "quench/ModelAux.h"
#include "quench/Model.h"
#include "quench/ModelTrajectory.h"
#include "quench/Geometry.h"
#include "quench/State.h"
#include "quench/Variables.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class StateTestFixture : TestFixture {
 public:
  StateTestFixture() {
    file_.reset(new eckit::LocalConfiguration(TestConfig::config(), "state"));
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    date_str_ = file_->getString("date");
    time_.reset(new util::DateTime(date_str_));
    vars_.reset(new quench::Variables(TestConfig::config()));
  }
  ~StateTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> file_;
  quench::Model * model_;  // Not actually used
  std::unique_ptr<quench::Geometry> resol_;
  std::string date_str_;
  std::unique_ptr<util::DateTime> time_;
  std::unique_ptr<quench::Variables> vars_;
};
// -----------------------------------------------------------------------------
CASE("test_State") {
  StateTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_State_readin_constructor") {
    std::unique_ptr<quench::State>
      xx(new quench::State(*fix.resol_, *fix.model_, *fix.file_));
    EXPECT(xx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_interpolation_constructor") {
    quench::State xx1(*fix.resol_, *fix.model_, *fix.file_);
    std::unique_ptr<quench::State> xx2(new quench::State(*fix.resol_, *fix.model_, xx1));
    EXPECT(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_copy_constructor") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);
    std::unique_ptr<quench::State> xx2(new quench::State(xx));
    EXPECT(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_destructor") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_getField") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);

    // there are 2 values in FieldL95: the 1st is the *resol_ value,
    // the 2nd is a vector of doubles initialised to 0.0, the size of the
    // vector is the *resol_ value
    EXPECT(xx.getField().resol() == fix.resol_->npoints());
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_validTime") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);
    EXPECT(xx.validTime().toString() == fix.date_str_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_assignment") {
    quench::State xx1(*fix.resol_, *fix.model_, *fix.file_);

    // construct the second State object
    quench::State xx2(*fix.resol_, *fix.model_, *fix.file_);

    xx2 = xx1;

    // initially, xx1 held data, xx2 was initialised to 0, now
    // ensure the two States are the same
    EXPECT(xx1.validTime() == xx2.validTime());
    EXPECT(xx1.getField().resol() == xx2.getField().resol());
    for (int i = 0; i < xx1.getField().resol(); ++i) {
      EXPECT(xx1.getField()[i] == xx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_compound_assignment") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);

    // construct the Increment object
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    xx += dx;

    // both xx and dx were initialised with the same data,
    // so when the two are added together xx should be double incL95
    for (int i = 0; i < xx.getField().resol(); ++i) {
      EXPECT(xx.getField()[i] == 2.0 * dx.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_read") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);

    // to verify the information has been read correctly, we need to open
    // and read the file using ifstream functionality
    const std::string filename(fix.file_->getString("filename"));
    std::ifstream inStream(filename.c_str());
    if (!inStream.is_open()) {
      oops::Log::error() << "read functionality cannot be determined" << std::endl;
    }

    int resolInt;
    inStream >> resolInt;

    std::string time;
    inStream >> time;

    std::vector<double> doubleVec(resolInt);
    for (int i = 0; i < resolInt; ++i) {
      inStream >> doubleVec[i];
    }
    inStream.close();

    for (int i = 0; i < fix.resol_->npoints(); ++i) {
      EXPECT((xx.getField())[i] == doubleVec[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_State_timeUpdate") {
    quench::State xx(*fix.resol_, *fix.model_, *fix.file_);

    util::Duration dur(10);
    util::DateTime updatedTime(*fix.time_);
    updatedTime += dur;
    std::string updated_date_string(updatedTime.toString());

    xx.validTime() += dur;

    EXPECT(xx.validTime().toString() == updated_date_string);
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_State_stream_output") {
    quench::State xx(*fix.resol_, *fix.vars_, *fix.time_);
    xx.read(*fix.file_);

    // use the operator<< method to write the value to a file

    std::filebuf fb;
    std::string filename("StateTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << xx;
    fb.close();

    std::vector<double> doubleVec(xx.getField().resol());
    for (int i = 0; i < xx.getField().resol(); ++i) {
      doubleVec[i] = xx.getField()[i];
    }
    std::vector<double>::iterator iter;
    iter = std::min_element(doubleVec.begin(), doubleVec.end());
    double min = *iter;
    iter = std::max_element(doubleVec.begin(), doubleVec.end());
    double max = *iter;
//    double avg = 0.0;

    // then read the values that were written to the file
    int idxStart;
    int idxEnd;
    double inputDouble = 0.0;
    std::string input;
    std::string inputTest;
    std::string inputTimeTest(" Valid time: " + fix.file_->getString("date"));
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);  // ignore the first (blank) line
      getline(inputFile, input);

      EXPECT(input == inputTimeTest);

      getline(inputFile, input);
      // test min value
      idxStart = input.find("=");
      idxEnd = input.find(",");
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = eckit::Translator<std::string, double>()(inputTest);
      EXPECT(oops::is_close(inputDouble, min, 0.0001));
      // test max value
      idxStart = input.find("=", idxEnd);
      idxEnd = input.find(",", idxEnd + 1);
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = eckit::Translator<std::string, double>()(inputTest);
      EXPECT(oops::is_close(inputDouble, max, 0.0001));
      // test avg value
//      idxStart = input.find("=", idxEnd + 1);
//      idxEnd = input.find(",", idxEnd + 1);
//      inputTest = input.substr(idxStart + 1);
//      inputDouble = eckit::Translator<std::string, double>()(inputTest);
//      EXPECT(oops::is_close(inputDouble, avg, 0.0001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
*/
// -----------------------------------------------------------------------------
//  SECTION("test_State_step_traj") {
//    quench::State xx(*fix.resol_, *fix.vars_, *fix.time_);
//    xx.read(*fix.file_);
//
//    // construct a Model object
//    eckit::LocalConfiguration modelCfg(TestConfig::config(), "model");
//    quench::Model Model(*fix.resol_, modelCfg);
//
//    // construct a ModelAux object
//    eckit::LocalConfiguration biasCfg(TestConfig::config(), "ModelAux");
//    quench::ModelAux ModelAux(*fix.resol_, biasCfg);
//
//    // construct a ModelTrajectory object
//    quench::ModelTrajectory modelTraj(true);
//    // copy construct a FieldL95 object to set params
//    // in the ModelTrajectory object
//    quench::FieldL95 fieldL95(xx.getField());
//    modelTraj.set(fieldL95);
//
//    xx.stepTraj(Model, ModelAux, modelTraj);
//
//    // create test data
//    Model.stepRK(fieldL95, ModelAux, modelTraj);
//
//    // test xx data against newly created fieldL95 test data
//    for(int i = 0; i < xx.getField().resol(); ++i) {
//      EXPECT(xx.getField()[i] == fieldL95[i]);
//    }
//  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
