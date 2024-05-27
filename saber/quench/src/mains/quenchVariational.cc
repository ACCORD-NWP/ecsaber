/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "quench/instantiateQuenchMatrices.h"
#include "quench/Logbook.h"
#include "quench/Traits.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  quench::instantiateQuenchMatrices();
  saber::instantiateCovarFactory<quench::Traits>();
  oops::Variational<quench::Traits> var;
  quench::Logbook::start();
  run.execute(var);
  quench::Logbook::stop();
  return 0;
}
