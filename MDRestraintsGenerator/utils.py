"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains utility functions.
"""

import MDAnalysis as mda
import warnings


def _get_host_anchors(atomgroup, l_atom, anchor_selection, num_atoms=3,
                     init_cutoff=5, max_cutoff=9):
    """Tool to get the host anchor carbon based nearest to l_atom.

    This function will first start by finding available atoms within the
    init_cutoff distance (default 5 Å). If less than num_atoms are found,
    it will expand by 1 Å until it hits max_cutoff.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    l_atom : int
        index of ligand atom to search around.
    anchor_selection : str
        Selection string for atoms to include in the search.
    num_atoms : int
        Minimum number of anchor atoms to return [3]
    init_cutoff : float
        Initial cutoff (Å) to search for anchor atoms in [5]
    max_cutoff : float
        Maximum cutoff (Å) to search for anchor atoms in [9]


    Notes
    -----
    Currently this only searches on the first frame, in the future this will
    likely be expanded to cover all frames.
    """

    found_atoms = 0
    cutoff = init_cutoff

    while found_atoms < num_atoms:
        sel_str = f"{anchor_selection} and around {cutoff} index {l_atom}"
        anchor_atoms = atomgroup.select_atoms(sel_str)
        found_atoms = anchor_atoms.n_atoms
        if found_atoms < num_atoms:
            if cutoff < max_cutoff:
                cutoff += 1
                wmsg = (f"too few anchor atoms found, expanding cutoff to "
                        f"{cutoff} Å")
                warnings.warn(wmsg)
            else:
                errmsg = "Too few anchor type atoms found"
                raise RuntimeError(errmsg)

    return [ix for ix in anchor_atoms.ix]


def _get_bonded_cn_atoms(atomgroup, anchor_ix):
    """Helper function to get the neighbouring atoms of a CA atom

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    anchor_ix : int
        Index of the CA atom residue for which to pick up C and N atoms
    """
    anchor_ag = atomgroup.select_atoms(f'index {anchor_ix}')
    resid = anchor_ag.resid
    c_atom_str = f"(backbone) and (name C) and (resid {resid})"
    n_atom_str = f"(backbone) and (name N) and (resid {resid+1})"

    c_atom = atomgroup.select_atoms(c_atom_str)
    n_atom = atomgroup.select_atoms(n_atom_str)

    return c_atom.atoms[0].ix, n_atom.atoms[0].ix


def _get_bonded_atoms(atomgroup, anchor_ix):
    """Helper function to get heavy atoms bonded to anchor

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    anchor_ix : int
        Index of the CA atom residue for which to pick up C and N atoms
    """
    # First get all heavy atoms bonded to the anchor (and not hydrogen)
    sel_str = f"(bonded index {anchor_ix}) and not name H*"
    atg1 = atomgroup.select_atoms(sel_str)
    
    # Loop over found indexes and see if any of them have imoqie bpmded atoms
    for index in [ix for ix in atg.ix]:
        sel_str = (f"(bonded index {ix}) and not (index {anchor_ix}) and not "
                   f"name H*")
        atg2 = atomgroup.select_atoms(sel_str)
        if atg2.n_atoms > 0:
            break

    try:
        return index, atg2.atoms[0].ix
    except IndexError:
        errmsg = "could not find p_atom[2]"
        raise RuntimeError(errmsg) from None


def get_host_atoms(atomgroup, l_atom, p_selection, num_restraints=3,
                   protein_routine=True, search_init_cutoff=5,
                   search_max_cutoff=9):
    """Gets a list of restraint host atoms using a distance search around
    the ligand bonding atoms.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    l_atom : int
        Index of the ligand atom involved in the bond formation.
    p_selection : str
        Selection string to define the atoms to be included as possible
        anchor points (direct binding to `l_atom`). Note: "protein and name CA"
        will default to a sepcial routine looking for C/N atoms bonded to CA
        unless `protein_routine` is `False`.
    num_restraints : int
        Number of boresch restraints to try to generate for a given ligand
        anchor atom. [3]
    protein_routine : bool
        Option to turn off the C/N atom gathering routine if `p_selection` is
        passed as "protein and name CA"
    search_init_cutoff : float
        Minimum cutoff distance to look for host anchor atoms. [5.0]
    search_max_cutoff : float
        Maximum cutoff distance to look for host anchor atoms. [9.0]
    """
    host_anchors = _get_host_anchors(atomgroup, l_atom, anchor_selection,
                                     num_atoms=num_restraints,
                                     init_cutoff=search_init_cutoff,
                                     max_cutoff=search_max_cutoff)

    p_atoms = [ [i] for i in host_anchors ]

    for entry in p_atoms:
        if protein_routine and (p_selection == "protein and name CA"):
            ix2, ix3 = _get_bonded_cn_atoms(atomgroup, entry[0])
        else:
            ix2, ix3 = _get_bonded_atoms(atomgroup, entry[0])

        entry.extend([ix2, ix3])

    return p_atoms

