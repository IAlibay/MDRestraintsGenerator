"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains IO functions.
"""


def _write_bond_header(rfile):
    """Helper function to write the bond intermolecular section"""
    rfile.write('[ bonds ]\n')
    rfile.write(';    ai   aj    type  bA         kA       bB    kB\n')


def _write_bond(bond, index, force_constant, rfile):
    """Helper function to write out the bond part of the boresch restraint

    Parameters
    ----------
    bond : Bond object
    index : int
        index of the frame to be used
    force_constant : float
        force constant for bond strength (kcal/mol)
    rfile : file object
        output restraint file
    """
    # Get index values
    atom1 = bond.atomgroup.atoms[0].ix + 1
    atom2 = bond.atomgroup.atoms[1].ix + 1
    length = bond.values[index] / 10.0
    bond_fc = force_constant * 4.184 * 100
    rfile.write(f"{atom1:>6}{atom2:>6}    6    {length:>6.3f}      0.0   "
                f"{length:>6.3f}   {bond_fc:>6.2f}\n")


def _write_angle_header(rfile):
    """Helper function to write the angle intermolecular section"""
    rfile.write('[ angles ]\n')
    rfile.write(';   ai    aj    ak    type    thA      '
                'fcA      thB       fcB\n')


def _write_angle(angle, index, force_constant, rfile):
    """Helper function to write out the angle part of the boresch restraint

    Parameters
    ----------
    angle : Angle object
    index : int
        index of the frame to be used
    force_constant : float
        force constant for the angle restraint
    rfile : file object
        output restraint file
    """
    angle_fc = force_constant * 4.184

    atom1 = angle.atomgroup.atoms[0].ix + 1
    atom2 = angle.atomgroup.atoms[1].ix + 1
    atom3 = angle.atomgroup.atoms[2].ix + 1
    val = angle.values[index]

    rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}       1   {val:>6.3f}"
                f"    0.0   {val:>6.3f}    {angle_fc:>6.2f}\n")


def _write_dihedral_header(rfile):
    """Helper function to write the dihedral intermolecular section"""
    rfile.write('[ dihedrals ]\n')
    rfile.write(';   ai    aj    ak    al  type    phiA      fcA    '
                'phiB      fcB\n')


def _write_dihedral(dihedral, index, force_constant, rfile):
    """Helper function to write out the dihedral part of the boresch restraint

    Parameters
    ----------
    dihedral : Dihedral object
    index : int
        index of the frame to be used
    force_constant : float
        force constant for the dihedral restraint
    rfile : file object
        output restraint file
    """
    angle_fc = force_constant * 4.184

    atom1 = dihedral.atomgroup.atoms[0].ix + 1
    atom2 = dihedral.atomgroup.atoms[1].ix + 1
    atom3 = dihedral.atomgroup.atoms[2].ix + 1
    atom4 = dihedral.atomgroup.atoms[3].ix + 1
    val = dihedral.values[index]
    rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}{atom4:>6}     2    "
                f"{val:>6.3f}    0.0    {val:>6.3f}   {angle_fc:>6.2f}\n")
