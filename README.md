# ECSABER
Interface between OOPS (ECMWF) and SABER (JEDI).

## JEDI sources processing
ECSABER downloads the code of three JEDI repositories as submodules:
- `oops-jedi`: https://github.com/JCSDA/oops
- `vader-jedi`: https://github.com/JCSDA/vader
- `saber-jedi`: https://github.com/JCSDA/saber

Some `oops-jedi` source files are copied and modified to create the `oops-patch` library, which contains missing functionalities of OOPS (ECMWF).

Similarly, the `vader-jedi` and `saber-jedi` source files are copied and modified to create the `saber` and `vader` libraries, respectively.

## Installation
To install ECSABER
1) Create a source directory `$SRC_DIR` and copy the following code into a file named `CMakeLists.txt`:
```
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.6 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project( oops-bundle VERSION 0.0.1 LANGUAGES C CXX Fortran )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

# Enable OpenMP and MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )

# Define bundle
ecbuild_bundle_initialize()

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# when building, already use the install RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# ECMWF tools
ecbuild_bundle( PROJECT eckit    GIT "https://github.com/ecmwf/eckit.git" TAG 1.24.4 )
ecbuild_bundle( PROJECT fckit    GIT "https://github.com/ecmwf/fckit.git" TAG 0.11.0 )
ecbuild_bundle( PROJECT fiat     GIT "https://github.com/ecmwf-ifs/fiat.git" BRANCH main )
ecbuild_bundle( PROJECT ectrans  GIT "https://github.com/ecmwf-ifs/ectrans.git" BRANCH main )
ecbuild_bundle( PROJECT atlas    GIT "https://github.com/ecmwf/atlas.git" TAG 0.35.0 )

# OOPS
ecbuild_bundle( PROJECT oops     GIT "https://git.ecmwf.int/scm/oops/oops.git" BRANCH develop UPDATE )

# ECSABER
ecbuild_bundle( PROJECT ecsaber  GIT "https://git.ecmwf.int/scm/sab/ecsaber.git" BRANCH develop UPDATE RECURSIVE )

ecbuild_bundle_finalize()
```
2) Create a build director `$BLD_DIR`.
3) Go to `$BLD_DIR`.
4) Configure: `ecbuild $SRC_DIR`
5) Compile: `make -jN` where N is the number of threads
6) Test: `ctest`

## Contact
Benjamin Menetrier, Meteorologisk Institutt, Norway.<br>
benjamin.menetrier -at- met.no
