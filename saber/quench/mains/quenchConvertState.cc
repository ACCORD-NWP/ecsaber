/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <limits>

#include "oops/runs/ConvertState.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"
#include "src/instantiateQuenchMatrices.h"
#include "src/Logbook.h"
#include "src/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::Log::test().setf(std::ios::scientific);
  oops::Log::test().precision(std::numeric_limits<double>::digits10+1);
  quench::instantiateQuenchMatrices();
  oops::ConvertState<quench::Traits> cs;
  quench::Logbook::start();
  run.execute(cs);
  quench::Logbook::stop();
  return 0;
}
