"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains IO functions.
"""

import MDAnalysis as mda
import warnings


def _write_boresch_bond(bond, index, force_constant, rfile):
    """Helper function to write out the bond part of the boresch restraint

    Parameters
    ----------
    bond : Bond object
    force_constant : float
        force constant for bond strength
    rfile : file object
        output restraint file
    """
    rfile.write('[ bonds ]\n')
    rfile.write(';    ai   aj    type  bA         kA       bB    kB\n')
    # Get index values
    atom1 = bond.atomgroup.atoms[0].ix + 1
    atom2 = bond.atomgroup.atoms[1].ix + 2
    length = bond.values[index] / 10.0
    bond_fc = force_constant * 4.184 * 100
    rfile.write(f"{atom1:>6}{atom2:>6}    6    {length:>6.3f}      0.0   "
                f"{length:>6.3f}   {bond_fc:>6.2f}\n")


def _write_boresch_angles(angles, index, force_constant, rfile):


def _write_boersch_dihedrals(dihedrals, index, force_constant, rfile):


def write_boresch(atomgroup, restraint, index, force_constant=10.0):
    """Tool to get the host anchor carbon based nearest to l_atom.

    This function will first start by finding available atoms within the
    init_cutoff distance (default 5 Å). If less than num_atoms are found,
    it will expand by 1 Å until it hits max_cutoff.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to create a restraint.
    restraint : BoreschRestraint
        Object containing the boresch restraint information.
    index : int
        Index of frame to be written out.
    force_constant : float
        Force constant to use for the restraint.

    Note
    ----
    For now only gmx BoreschRestraints are supported, in the future
    we hope to change this functionality to allow for various MD engines.
    """
    bond = restraint.bond.values[index]
    angle1 = restraint.angles[0].values[index]
    angle2 = restraint.angles[1].values[index]
    dih1 = restraint.dihedrals[0].values[index]
    dih2 = restraint.dihedrals[1].values[index]
    dih3 = restraint.dihedrals[2].values[index]

    # Generate gro file
    ag = atomgroup.atoms # in case it's a universe
    ag.universe.trajectory[index]
    ag.write('ClosestRestraintFrame.gro')

    with open('BoreschRestraint.top', 'w') as rfile:
        rfile.write('\n')
        rfile.write('; restraints\n')
        rfile.write('[ intermolecular_interactions ]\n')

        _write_boresch_bond(restraint.bond, index, force_constant, rfile)
        _write_boresch_angles(restraint.angles, index, force_constant, rfile)
        _write_boresch_dihedrals(restraint.dihedrals, index, force_constant,
                                rfile)
