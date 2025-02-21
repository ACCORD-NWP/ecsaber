/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_
#define OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_

#include "oops/assimilation/GETKFSolver.h"
//#include "oops/assimilation/GETKFSolverPert.h"
#include "oops/assimilation/LETKFSolver.h"
//#include "oops/assimilation/LETKFSolverGSI.h"
//#include "oops/assimilation/LETKFSolverPert.h"
#include "oops/assimilation/LocalEnsembleSolver.h"

namespace oops {

template <typename MODEL> void instantiateLocalEnsembleSolverFactory() {
  static LocalEnsembleSolverMaker<MODEL, LETKFSolver<MODEL> > makerLETKF_("LETKF");
//  static LocalEnsembleSolverMaker<MODEL, LETKFSolverGSI<MODEL> > makerGSI_("GSI LETKF");
//  static LocalEnsembleSolverMaker<MODEL, LETKFSolverPert<MODEL>
//         > makerLETKFPert_("Perturbed LETKF");
  static LocalEnsembleSolverMaker<MODEL, GETKFSolver<MODEL> > makerGETKF_("GETKF");
//  static LocalEnsembleSolverMaker<MODEL, GETKFSolverPert<MODEL>
//         > makerGETKFPert_("Perturbed GETKF");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_
