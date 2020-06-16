MDRestraintsGenerator
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/MDRestraintsGenerator.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/MDRestraintsGenerator)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/MDRestraintsGenerator/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/MDRestraintsGenerator/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/MDRestraintsGenerator/branch/master)

# MDRestraintsGenerator

A framework for generating restraints for MD simulations (from MD simulations).

The code currently only looks at BoreschRestraints, but we aim to extend to:

- Hardwall restraints
- Harmonic restraints
- Attach Pull Restraint style restraints
- Complex multidimensional restraints

Note: This is considered to be pre-production code, a lot is left to do and
major changes will happen at any time.

## How to use

Currently implemented; FindBoreschRestraint

This works as a usual AnalysisBase object, please see `scripts` for a basic
usage example.

### Dependencies

- MDAnalysis
- NumPy
- SciPy
- Matplotlib

### Copyright

Copyright (c) 2020, Irfan Alibay


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
