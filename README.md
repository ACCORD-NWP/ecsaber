# ECSABER
Interface between OOPS (ECMWF) and SABER (JEDI).

## Installation
To install ECSABER
1) Clone this code into a source directory `$SRC_DIR`.
2) Create a build director `$BLD_DIR`.
3) Go to `$BLD_DIR`.
4) Configure: `ecbuild $SRC_DIR/bundle`
5) (optional) Go to the `ecsaber`
6) Compile: `make -jN` where N is the number of threads
7) Test: `ctest`

## JEDI sources processing
If the `ECSABER_UPDATE` environment variable is set to `ON`, then ECSABER will update its source code from three JEDI repositories:
- `oops-jedi`: https://github.com/JCSDA/oops
- `saber-jedi`: https://github.com/JCSDA/saber
- `vader-jedi`: https://github.com/JCSDA/vader

The `oops-jedi`, `saber-jedi` and `vader-jedi` repository should be cloned manually in $SRC_DIR/ecsaber, or linked to other source directories. Then at configure time, `ecbuild` will prompt the user with source update options.

## Contact
Benjamin Menetrier, Meteorologisk Institutt, Norway.<br>
benjamin.menetrier -at- met.no
