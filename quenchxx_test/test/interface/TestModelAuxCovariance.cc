/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "oops/runs/Run.h"
#include "quench/Traits.h"
#include "test/interface/ModelAuxCovariance.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ModelAuxCovariance<quench::Traits> tests;
  run.execute(tests);
  return 0;
}
