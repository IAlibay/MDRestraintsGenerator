"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains IO functions.


.. versionchanged:: 0.2.0
   [intermolecular_interactions] writer helper functions have been renamed with
   the `intinters` tag to denote their use.

"""


def yes_no(option):
    """Helper function to return yes/no based on a bool entry


    .. versionadded:: 0.2.0
    """
    return 'yes' if option else 'no'


def _write_intinters_header(rfile):
    """Helper function to write the header for an intermolecular_interactions
    section of a GMX topology file


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    rfile.write('\n')
    rfile.write('; restraints\n')
    rfile.write('[ intermolecular_interactions ]\n')


def _write_intinters_bond_header(rfile):
    """Helper function to write the bond intermolecular section


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    rfile.write('[ bonds ]\n')
    rfile.write(';    ai   aj    type  bA         kA       bB    kB\n')


def _write_intinters_bond(bond, index, force_constant, rfile):
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



    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    # Get index values
    atom1 = bond.atomgroup.atoms[0].ix + 1
    atom2 = bond.atomgroup.atoms[1].ix + 1
    length = bond.values[index] / 10.0
    bond_fc = force_constant * 4.184 * 100
    rfile.write(f"{atom1:>6}{atom2:>6}    6    {length:>6.3f}      0.0   "
                f"{length:>6.3f}   {bond_fc:>6.2f}\n")


def _write_intinters_angle_header(rfile):
    """Helper function to write the angle intermolecular section


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    rfile.write('[ angles ]\n')
    rfile.write(';   ai    aj    ak    type    thA      '
                'fcA      thB       fcB\n')


def _write_intinters_angle(angle, index, force_constant, rfile):
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


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    angle_fc = force_constant * 4.184

    atom1 = angle.atomgroup.atoms[0].ix + 1
    atom2 = angle.atomgroup.atoms[1].ix + 1
    atom3 = angle.atomgroup.atoms[2].ix + 1
    val = angle.values[index]

    rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}       1   {val:>6.3f}"
                f"    0.0   {val:>6.3f}    {angle_fc:>6.2f}\n")


def _write_intinters_dihedral_header(rfile):
    """Helper function to write the dihedral intermolecular section


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    rfile.write('[ dihedrals ]\n')
    rfile.write(';   ai    aj    ak    al  type    phiA      fcA    '
                'phiB      fcB\n')


def _write_intinters_dihedral(dihedral, index, force_constant, rfile):
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


    .. versionadded:: 0.1.0
    .. versionchanged:: 0.2.0
       Renamed to include `intinters` in the name.
    """
    angle_fc = force_constant * 4.184

    atom1 = dihedral.atomgroup.atoms[0].ix + 1
    atom2 = dihedral.atomgroup.atoms[1].ix + 1
    atom3 = dihedral.atomgroup.atoms[2].ix + 1
    atom4 = dihedral.atomgroup.atoms[3].ix + 1
    val = dihedral.values[index]
    rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}{atom4:>6}     2    "
                f"{val:>6.3f}    0.0    {val:>6.3f}   {angle_fc:>6.2f}\n")


def _write_pull_header(rfile, pull=True, pull_print_com=False,
                       pull_print_ref_value=False, pull_print_components=False,
                       pull_nstxout=50, pull_nstfout=0,
                       pull_pbc_ref_prev_step_com=True,
                       pull_xout_average=False, pull_fout_average=False):
    """Helper function to write the header portion of an MDP pull section

    Disclaimer: MDP options are detailed as written in the GROMACS manual
    version 2021. See here for more details:
    https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html

    Parameters
    ----------
    rfile : file object
        output pull restraint MDP file
    pull : bool [`True`]
        Enable COM pulling (if `False`, all pull MDP options will be ignored)
    pull_print_com : bool [`False`]
        Turns on the printing of COM of all groups for all pull coordinates
    pull_print_ref_value : bool [`False`]
        Print the reference value for each pull coordinate
    pull_print_components : bool [`False`]
        Print the distance and Cartesian components selection in pull-coord-dim
    pull_nstxout : int [0]
        Frequency for writing out COMs for all the pull groups (0 is never)
    pull_nstfout : int [0]
        Frequency for writing out the force of all the pulled group
        (0 is never)
    pull_pbc_ref_prev_step_com : bool [`True`]
        Use the COM of the previous step as the reference for the treatment of
        periodic boundary conditions.
    pull_xout_average : bool [`False`]
        Write average coordinates for all the pull groups instead of
        instantaneous coordinates.
    pull_fout_average : bool [`False`]
        Write the average force instead of the instantaneous one for all the
        pulled groups.


    .. versionadded:: 0.2.0
    """

    # header
    rfile.write(';----------------------------------------------------\n')
    rfile.write('; PULL RESTRAINT OPTIONS\n')
    rfile.write(';----------------------------------------------------\n')

    # go through the options
    rfile.write(f"pull = {yes_no(pull)}\n")
    rfile.write(f"pull-print-com = {yes_no(pull_print_com)}\n")
    rfile.write(f"pull-print-ref-value = {yes_no(pull_print_ref_value)}\n")
    rfile.write(f"pull-print-components = {yes_no(pull_print_components)}\n")
    rfile.write(f"pull-nstxout = {pull_nstxout}\n")
    rfile.write(f"pull-nstfout = {pull_nstfout}\n")
    rfile.write(
        f"pull-pbc-ref-prev-step-com = {yes_no(pull_pbc_ref_prev_step_com)}\n")
    rfile.write(f"pull-xout-average = {yes_no(pull_xout_average)}\n")
    rfile.write(f"pull-fout-average = {yes_no(pull_fout_average)}\n")


def _write_pull_groups(rfile, group_names, group_pbc_atoms):
    """Helper function to write the group portion of an MDP pull section

    Parameters
    ----------
    rfile : file object
        output pull restraint MDP file
    group_names : list
        List of names for all pull groups.
    group_pbc_atoms : list
        Reference atoms used for the treatment of PBC conditions.

    Raises
    ------
    ValueError
        If group_names and group_pbc_atoms are of different length.

    TODO
    ----
    * Currently does not handle `pull-group?-weights`. If these are useful they
      can be implemented at a later stage.


    .. versionadded:: 0.2.0
    """

    if len(group_names) != len(group_pbc_atoms):
        errmsg = 'group names has a different length than group pbc atoms'
        raise ValueError(errmsg)

    # number of groups
    rfile.write(f"pull-ngroups = {len(group_names)}\n")

    # loop through names, index is 1-indexed
    for index, name in enumerate(group_names):
        rfile.write(f"pull-group{index+1}-name = {name}\n")

    # loop through pbc atoms, index and atom is 1-indexed
    for index, atom in enumerate(group_pbc_atoms):
        rfile.write(f"pull-group{index+1}-pbcatom = {atom+1}\n")


def _write_pull_coords(rfile,
                       pull_coord_types,
                       pull_coord_geometries,
                       pull_coord_groups,
                       pull_coord_dims,
                       pull_coord_starts,
                       pull_coord_inits,
                       pull_coord_rates,
                       pull_coord_ks,
                       pull_coord_kBs,
                       pull_coord_potential_providers=None,
                       pull_coord_origins=None,
                       pull_coord_vecs=None):
    """Helper function to write the coord portion of an MDP pull section

    Parameters
    ----------
    pull_coord_types : list of strings
        List of pull coord types for each pull coordinates.
    pull_coord_geometries : list of strings
        List of types of COM restraint goemetries.
    pull_coord_groups : list of tuples
        List of tuples containing the group indices (1-indexed) involved in
        each pull coordinate.
    pull_coord_dims : list of strings
        List of strings containing the cartesian coordinates on which the
        pull coordinate will act. Note: expected to change to a list of of
        bools in the future.
    pull_coord_starts : list of bools
        List of bools detailing whether or not to add the COM of the starting
        conformation to pull-coordN-init.
    pull_coord_inits : list of floats
        List of floats detailing the reference distance/angle at the start of
        the simulation, units: `nm` or `degrees`, formatted to %.3f.
    pull_coord_rates : list of floats
        List of floats detailing the rate of change of the reference distance/
        angle in nm/ps or degrees/ps, Formatted to %.3f.
    pull_coord_ks : list of floats
        List of floats containing the force constants of each pull coordinate.
        Units in kJ mol^-1 nm^-2, kJ mol^-1 nm^-1, kJ mol^-1 rad^-2, or
        kJ mol^-1 rad^-1. Formatted to %.3f.
    pull_coords_kBs : list of floats
        List of floats as pull_coords_ks, but for state B when doing a free
        energy perturbation. Formatted to %.3f.
    pull_coord_potential_providers : list of strings [`None`]
        List of external potential providers, default `None` does not set this
        MDP value.
    pull_coord_origins : list of strings [`None`]
        Containing the pull reference position (format 'float, float, float')
    pull_coord_vecs : list of strings [`None`]
        Containing the pull direction (format 'float, float, float')

    TODO
    ----
    * Change pull_coord_dims, pull_coord_origins, and pull_coord_strings, to
      be using bools and floats instead of strings.


    .. versionadded:: 0.2.0
    """

    # number of coords
    rfile.write(f"pull-ncoords = {len(pull_coord_types)}\n")

    for i in range(len(pull_coord_types)):
        index = i + 1  # index is 1 formatted
        rfile.write(f"pull-coord{index}-type = {pull_coord_types[i]}\n")
        rfile.write(
            f"pull-coord{index}-geometry = {pull_coord_geometries[i]}\n")
        groups_string = ' '.join(map(str, pull_coord_groups[i]))
        rfile.write(
            f"pull-coord{index}-groups = {groups_string}\n")
        rfile.write(f"pull-coord{index}-dim = {pull_coord_dims[i]}\n")
        rfile.write(
            f"pull-coord{index}-start = {yes_no(pull_coord_starts[i])}\n")
        rfile.write(f"pull-coord{index}-init = {pull_coord_inits[i]:.3f}\n")
        rfile.write(f"pull-coord{index}-rate = {pull_coord_rates[i]:.3f}\n")
        rfile.write(f"pull-coord{index}-k = {pull_coord_ks[i]:.3f}\n")
        rfile.write(f"pull-coord{index}-kB = {pull_coord_kBs[i]:.3f}\n")

        # optional entries
        if pull_coord_potential_providers is not None:
            entry = (f"pull-coord{index}-potential-provider = "
                     f"{pull_coord_potential_providers[i]}\n")
            rfile.write(entry)

        if pull_coord_origins is not None:
            entry = (f"pull-coord{index}-origin = {pull_coord_origins[i]}\n")
            rfile.write(entry)

        if pull_coord_vecs is not None:
            entry = (f"pull-coord{index}-vec = {pull_coord_vecs[i]}\n")
            rfile.write(entry)
