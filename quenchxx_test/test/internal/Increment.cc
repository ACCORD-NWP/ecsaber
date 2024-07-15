/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <fstream>
#include <iostream>


#include "./TestConfig.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/runs/Test.h"

#include "util/DateTime.h"
#include "util/Logger.h"

#include "quench/GeoVaLs.h"
#include "quench/Increment.h"
#include "quench/Locations.h"
#include "quench/Variables.h"
#include "quench/Geometry.h"
#include "quench/State.h"

#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class IncrementTestFixture : TestFixture {
 public:
  IncrementTestFixture() {
    file_.reset(new eckit::LocalConfiguration(TestConfig::config(), "state"));
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new quench::Geometry(res));
    date_str_ = "2010-01-01T09:00:00Z";
    time_.reset(new util::DateTime(date_str_));
    vars_.reset(new quench::Variables(TestConfig::config()));
  }
  ~IncrementTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> file_;
  std::unique_ptr<quench::Geometry> resol_;
  quench::Model * model_; // Not actually used
  std::string date_str_;
  std::unique_ptr<util::DateTime> time_;
  std::unique_ptr<quench::Variables> vars_;
};
// -----------------------------------------------------------------------------
CASE("test_Increment") {
  IncrementTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_Increment_constructor") {
    std::unique_ptr<quench::Increment>
      dx(new quench::Increment(*fix.resol_, *fix.vars_, *fix.time_));
    EXPECT(dx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_interpolation_constructor") {
    std::unique_ptr<quench::Increment>
      dx1(new quench::Increment(*fix.resol_, *fix.vars_, *fix.time_));
    std::unique_ptr<quench::Increment>
      dx2(new quench::Increment(*fix.resol_, *dx1));
    EXPECT(dx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_copy_constructor") {
    std::unique_ptr<quench::Increment>
      dx1(new quench::Increment(*fix.resol_, *fix.vars_, *fix.time_));
    std::unique_ptr<quench::Increment> dx2(new quench::Increment(*dx1));
    EXPECT(dx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_diff") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    // construct the first State object
    quench::State xx1(*fix.resol_, *fix.model_, *fix.file_);

    // read in the state config info
    xx1.read(*fix.file_);
    quench::State xx2(xx1);

    // to vary the results a little, change the second State field values
    double fact = 0.75;
    dx *= fact;
    xx2 += dx;

    dx.diff(xx1, xx2);

    // to verify the diff method has worked correctly, we need
    // to open and read the file containing the Fields
    std::string filename(fix.file_->getString("filename"));
    std::ifstream inStream(filename.c_str());
    if (!inStream.is_open()) {
      oops::Log::error() << "diff functionality cannot be determined" << std::endl;
    }

    // we read in these two values but do not use them
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
      EXPECT(oops::is_close((dx.getField())[i],
                         doubleVec[i] - (doubleVec[i] + (doubleVec[i] * fact)),
                         1.0e-6));
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_zero") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    // first check that we have good data, ie, at least one element is non-zero
    bool goodData = false;
    for (int i = 0; i < fix.resol_->npoints() && goodData == false; ++i) {
      if ((dx.getField())[i] != 0) {
        goodData = true;
      }
    }

    if (!goodData) {
      oops::Log::error() <<
        "unable to test zero method, since test data is already all zero" << std::endl;;
    } else {
      dx.zero();

      for (int i = 0; i < fix.resol_->npoints(); ++i) {
        EXPECT(dx.getField()[i] == 0);
      }
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_zero_set_datetime") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    // first check that we have good data, ie, at least one element is non-zero
    bool goodData = false;
    for (int i = 0; i < fix.resol_->npoints() && goodData == false; ++i) {
      if (dx.getField()[i] != 0) {
        goodData = true;
      }
    }

    if (!goodData) {
      oops::Log::error() <<
        "unable to test zero method, since test data is already all zero" << std::endl;
    } else {
      const std::string modified_date_string("2010-01-01T10:35:00Z");
      const util::DateTime dtModified(modified_date_string);
      dx.zero(dtModified);

      for (int i = 0; i < fix.resol_->npoints(); ++i) {
        EXPECT(dx.getField()[i] == 0);
      }

      EXPECT(dx.validTime().toString() != fix.date_str_);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_assignment") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // construct the second dx object
    quench::Increment dx2(*fix.resol_, *fix.vars_, *fix.time_);
    double fact = 0.75;
    dx2.read(*fix.file_);
    dx2 *= fact;

    dx1 = dx2;

    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_compound_assignment_add") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // copy construct the second State object
    quench::Increment dx2(dx1);

    dx1 += dx2;

    // since the two Increment objects started off with the same data,
    // once they've been added together incL591 will be double what incL952 is
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == 2.0 * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_compound_assignment_subtract") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // copy construct the second State object
    quench::Increment dx2(dx1);

    dx1 -= dx2;

    // since the two Increment objects started off with the same data,
    // once incL952 has been subtracted from incL951, the result is zero
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_compound_assignment_multiply") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    // create a copy of the original data for testing against
    std::vector<double> testData(dx.getField().resol());
    for (unsigned int ii = 0; ii < testData.size(); ++ii) {
      testData.at(ii) = dx.getField()[ii];
    }

    double fact = 0.75;
    dx *= fact;

    for (int ii = 0; ii < dx.getField().resol(); ++ii) {
      EXPECT(dx.getField()[ii] == testData.at(ii) * fact);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_axpy") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // copy construct the second State object
    quench::Increment dx2(dx1);

    double fact = 0.75;
    dx1.axpy(fact, dx2);

    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i] + fact * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_dot_product_with") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // copy construct the second State object
    quench::Increment dx2(dx1);

    double dpwResult = dx1.dot_product_with(dx2);

    // prepare a value to test against
    double testResult = 0.0;
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      testResult += (dx1.getField()[i] * dx2.getField()[i]);
    }

    EXPECT(dpwResult == testResult);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_schur_product_with") {
    quench::Increment dx1(*fix.resol_, *fix.vars_, *fix.time_);
    dx1.read(*fix.file_);

    // copy construct the second State object
    quench::Increment dx2(dx1);

    dx1.schur_product_with(dx2);

    // both incL951 and incL952 started off with the same data in x_,
    // so to test incL951 against incL952xincL952 is a valid test
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i] * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  /*
  SECTION("test_Increment_interpolateTL") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);

    int maxVecSize = 5;

    // construct the Locations object
    std::vector<double> doubleVec(maxVecSize);
    // populate the vector with values from 0.1 to 0.5
    for (int i = 0; i < maxVecSize; ++i) {
      doubleVec[i] = ((i + 1) * 0.1);
    }
    quench::Locations Locations(doubleVec);

    // construct the GeoVaLs object
    std::vector<int> intVec(maxVecSize);
    quench::GeoVaLs GeoVaLs(intVec);
    // populate the locval_ vector with values from 1.1 to 0.5
    for (int i = 0; i < intVec.size(); ++i) {
      GeoVaLs[i] = ((i + 1) * 0.1);
    }

    int origCurrent = GeoVaLs.current();

    dx.interpolateTL(Locations, GeoVaLs);

    for (int i = 0; i < Locations.nobs(); ++i) {
      EXPECT(GeoVaLs[origCurrent + i] == dx.getField()[i]);
    }
    EXPECT(GeoVaLs.current() == (origCurrent + Locations.nobs()));
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_interpolateAD") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);

    int maxVecSize = 5;

    // construct the Locations object
    std::vector<double> doubleVec(maxVecSize);
    // populate the vector with values from 0.1 to 0.5
    for (int i = 0; i < maxVecSize; ++i) {
      doubleVec[i] = ((i + 1) * 0.1);
    }
    quench::Locations Locations(doubleVec);

    // construct the GeoVaLs object
    std::vector<int> intVec(maxVecSize);
    quench::GeoVaLs GeoVaLs(intVec);
    // populate the locval_ vector with values from 0.1 to 0.5
    for (int i = 0; i < intVec.size(); ++i) {
      GeoVaLs[i] = ((i + 1) * 0.1);
    }

    int origCurrent = GeoVaLs.current();

    dx.interpolateAD(Locations, GeoVaLs);

    const double dres = static_cast<double>(dx.getField().resol());
    if (GeoVaLs.current() == 0) {
      GeoVaLs.current() = GeoVaLs.nobs();
    }
    GeoVaLs.current() -= Locations.nobs();

    for (int i = 0; i < Locations->nobs(); ++i) {
      int idx = round(Locations[i] * dres);
      EXPECT(dx.getField()[idx] == GeoVaLs[GeoVaLs.current() + i]);
    }
  }
  */
// -----------------------------------------------------------------------------
  SECTION("test_Increment_read") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

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
      EXPECT(dx.getField()[i] == doubleVec[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_write") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    dx.read(*fix.file_);

    eckit::LocalConfiguration opFileCfg(TestConfig::config(), "outputFile");
    dx.write(opFileCfg);

    // Should read back in and compare values
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_validTime") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);
    EXPECT(dx.validTime().toString() == fix.date_str_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_Increment_timeUpdate") {
    util::Duration dur("PT1H30M10S");
    std::string updated_date_string("2010-01-01T10:30:10Z");

    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);

    dx.validTime() += dur;

    EXPECT(dx.validTime().toString() == updated_date_string);
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_Increment_stream_output") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("IncrementTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dx;
    fb.close();

    // then read the value that was written to the file
    std::string input;
    std::string inputTest(" Valid time: " + fix.date_str_);
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);  // ignore the first (blank) line
      getline(inputFile, input);

      EXPECT(input == inputTest);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
*/
// -----------------------------------------------------------------------------
  SECTION("test_Increment_getField") {
    quench::Increment dx(*fix.resol_, *fix.vars_, *fix.time_);

    // there are 2 values in FieldL95: the 1st is the *resol_ value,
    // the 2nd is a vector of doubles initialised to 0.0, the size of the
    // vector is the *resol_ value (we're just checking the final one)
    EXPECT(dx.getField().resol() == fix.resol_->npoints());
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
