MDRestraintsGenerator
==============================
[//]: # (Badges)
![GH Actions CI](https://github.com/IAlibay/MDRestraintsGenerator/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/IAlibay/MDRestraintsGenerator/branch/master/graph/badge.svg)](https://codecov.io/gh/Bigginlab/MDRestraintsGenerator/branch/master)
[![DOI](https://zenodo.org/badge/185426662.svg)](https://zenodo.org/badge/latestdoi/185426662)

# MDRestraintsGenerator

A framework for generating restraints for MD simulations (from MD simulations).

The code currently only looks at BoreschRestraints, but we aim to extend to:

- Hardwall restraints
- Harmonic restraints
- Attach Pull Restraint style restraints
- Complex multidimensional restraints

Note: This is non-mature code, a lot is left to do and
major changes will happen at any time.

## Installation

Installation is currently only possible from source. This can be done in the following manner:

```
git clone https://github.com/bigginlab/MDRestraintsGenerator.git
cd MDRestraintsGenerator
pip install .
```

## How to use

The code currently only implements a means of deriving Boresch restraints for GROMACS simulations.
To achieve this, the following underlying methods are provided:

  1) A function to pick stable points in ligands for restraint attachment
     (`search.find_ligand_atoms`).
  2) A class for picking host restraint addition points (`search.FindHostAtoms`).
  3) A class for analysing a list of possible Boresch restraints over a given MD simulation and
     finding the most stable choice of restraint atoms (`restraints.FindBoreschRestraint`).
     
Boresch restraints are implemented under the BoreschRestraint class. When using
`restraints.FindBoreschRestraint`, once run (using the `run()` method), the preffered restraint
will be stored as such an object under the `restraint` attribute. The BoreschRestraint class
offers three useful methods:
  1) The `plot()` function which outputs images of the distributions for the each component
     of the Boresch restraint (one bond, two angles, three dihedrals). In addition to the
     histograms, indicating both the mean and the picked frame positions, Q Q plots are
     also given to show how close to normality the distribution is. The latter can be useful
     when trying to work out if the chosen variable may occupy different binding orientations.
  2) The `write()` function, which writes out the `ClosestRestraint.gro` and
     `BoreschRestraint.top` files. These are based on the "picked frame", either user supplied
     or, in most cases, automatically obtained via the `FindBoreschRestraint` routine as the
     "frame closest to the mean across all bond/angle/dihedral distributions". This `.gro`
     file outputs the system at that frame, and the `.top` file contains the
     "intermoecular_interactions" section of a GROMACS `.top` file. This can then be pasted
     into an existing `.top` file to apply the restraint.
  3) The `standard_state` function, which currently resturns the analytical standard state
     correction for the restraint addition.

An example use script is provided under `scripts.BoreschRestraintGMX.py`. Documentation docstrings
are provided for all functions and classes. These can be accessed by calling `help(function)`.

## Testing

A set of unit tests are provided under `MDRestraintsGenerator.tests`. To run these you will need
to have `pytest` installed. The tests can be run in the following manner:
```
pytest -v MDRestraintsGenerator.tests
```

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
