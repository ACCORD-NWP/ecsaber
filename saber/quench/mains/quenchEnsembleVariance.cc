/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "saber/oops/EnsembleVariance.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "src/instantiateQuenchMatrices.h"
#include "src/Logbook.h"
#include "src/Traits.h"

int main(int argc, char** argv) {
  oops::Run run(argc, argv);
  quench::instantiateQuenchMatrices();
  saber::instantiateCovarFactory<quench::Traits>();
  saber::EnsembleVariance<quench::Traits> ev;
  quench::Logbook::start();
  run.execute(ev);
  quench::Logbook::stop();
  return 0;
}
