/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * s
 */

#include "oops/runs/Run.h"
#include "quench/Traits.h"
#include "test/interface/Locations.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::Locations<quench::Traits> tests;
  run.execute(tests);
  return 0;
}
