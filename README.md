MDRestraintsGenerator
==============================
[//]: # (Badges)
![GH Actions CI](https://github.com/IAlibay/MDRestraintsGenerator/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/IAlibay/MDRestraintsGenerator/branch/master/graph/badge.svg)](https://codecov.io/gh/Bigginlab/MDRestraintsGenerator/branch/master)
[![DOI](https://zenodo.org/badge/185426662.svg)](https://zenodo.org/badge/latestdoi/185426662)

# MDRestraintsGenerator

A framework for generating restraints for MD simulations (from MD simulations).

The code currently implements a means of deriving Boresch-style restraints,
with exporters for GROMACS. There is also experimental code for COM based
restraints (i.e. harmonic distance or hard wall), which export to the gromacs
pull code. These experimental implementations have yet to be completely tested.

In future implementations we aim to expand to other MD engines (notably OpenMM
support will be coming in the near future as part of efforts to support
the work done by OpenFE).

We also aim to eventually implement the following restraint types:

- Attach Pull Restraint style restraints
- Arbitrary multidimensional restraints (will require API overhaul)

Note: This is non-mature code, a lot is left to do and
major changes will happen at any time.

## Installation

Installation can either be done via PyPi or from source.

To install the latest release via PyPi do:

```
pip install MDRestraintsGenerator
```

Installing the latest development code from source can be done using the
following:

```
git clone https://github.com/bigginlab/MDRestraintsGenerator.git
cd MDRestraintsGenerator
pip install .
```

## How to use

The code currently focuses on implementing a means of deriving Boresch restraints for GROMACS simulations.
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

To cite this code, please add refer the following:

  - https://doi.org/10.26434/chemrxiv-2022-cw2kq-v3
  - https://doi.org/10.5281/zenodo.6972482
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
