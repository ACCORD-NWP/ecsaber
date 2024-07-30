/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "saber/oops/EvilUpdate.h"
#include "src/instantiateQuenchMatrices.h"
#include "src/Logbook.h"
#include "src/Traits.h"

int main(int argc, char** argv) {
  oops::Run run(argc, argv);
  quench::instantiateQuenchMatrices();
  saber::EvilUpdate<quench::Traits> eu;
  quench::Logbook::start();
  run.execute(eu);
  quench::Logbook::stop();
  return 0;
}
