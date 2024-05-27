/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/MakeObs.h"
#include "oops/runs/Run.h"
#include "quench/instantiateQuenchMatrices.h"
#include "quench/Logbook.h"
#include "quench/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  quench::instantiateQuenchMatrices();
  oops::MakeObs<quench::Traits> mo;
  quench::Logbook::start();
  run.execute(mo);
  quench::Logbook::stop();
  return 0;
}
