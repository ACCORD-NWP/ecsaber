/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/QgTraits.h"
#include "model/LogbookQG.h"
#include "oops/runs/EnsembleVariance.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::EnsembleVariance<qg::QgTraits> ev;
  qg::LogbookQG::start();
  run.execute(ev);
  qg::LogbookQG::stop();
  return 0;
};
