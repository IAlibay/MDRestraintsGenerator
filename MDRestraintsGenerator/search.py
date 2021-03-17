"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains utility functions.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.lib.distances import capped_distance
import warnings


def _search_from_capped(ref_ag, config_ag, universe, cutoff,
                        return_distances=True):
    """Helper function to search for neighbouring atoms using
    capped_distances

    Parameters
    ----------
    ref_ag: MDAnalysis.AtomGroup
        Reference atoms to search around.
    config_ag: MDAnalysis.AtomGroup
        Atoms to include in the search.
    universe: MDAnalysis.Universe
        Universe of atomgroups passed (used for box dimensions).
    cutoff: float
        Cutoff search distance value from reference.
    return_distance: bool [`True`]
        Toggles the return of distances for each atom pair.

    Returns
    -------
    results: list
        List of up to two items containing the atom pairs within the cutoff,
        and if return_distances is True, the distances for each pair.


    .. versionchanged:: 0.2.0
       Now optionally returns the distances.
       Reference AtomGroup can be more than one atom.
    """
    results = capped_distance(ref_ag.atoms.positions,
                              config_ag.atoms.positions,
                              max_cutoff=cutoff,
                              box=universe.dimensions,
                              return_distances=return_distances)

    return results


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


def _get_bonded_host_cn_atoms(atomgroup, anchor_ix):
    """Helper function to get the neighbouring atoms of a CA atom

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    anchor_ix : int
        Index of the CA atom residue for which to pick up C and N atoms
    """
    anchor_ag = atomgroup.select_atoms(f'index {anchor_ix}')
    resid = anchor_ag.atoms[0].resid
    c_atom_str = f"(backbone) and (name C) and (resid {resid})"
    c_atom = atomgroup.select_atoms(c_atom_str)
    c_atom_ix = c_atom.atoms[0].ix
    n_atom_str = (f"(backbone) and (name N) and (resid {resid+1}) and "
                  f"(bonded index {c_atom_ix})")
    n_atom = atomgroup.select_atoms(n_atom_str)

    if n_atom:
        return c_atom.atoms[0].ix, n_atom.atoms[0].ix
    else:
        # Issue 30 - terminal residue
        # reverse selection to go to residue-1 rather than residue+1
        n_atom_str = f"(backbone) and (name N) and (resid {resid})"
        n_atom = atomgroup.select_atoms(n_atom_str)
        n_atom_ix = n_atom.atoms[0].ix
        c_atom_str = (f"(backbone) and (name C) and (resid {resid-1}) and "
                      f"(bonded index {n_atom_ix})")
        c_atom = atomgroup.select_atoms(c_atom_str)
        return n_atom.atoms[0].ix, c_atom.atoms[0].ix


def _get_bonded_host_atoms(atomgroup, anchor_ix,
                           exclusion_str="and not name H*"):
    """Helper function to get heavy atoms bonded to anchor

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    anchor_ix : int
        Index of the CA atom residue for which to pick up bonded atoms
    exclusion_str : str
        Exclusion selection string to not pick up certain types of atoms.
        ['and not name H*']
    """
    # First get all heavy atoms bonded to the anchor (and not hydrogen)
    sel_str = f"(bonded index {anchor_ix}) {exclusion_str}"
    atg1 = atomgroup.select_atoms(sel_str)

    if len(atg1.atoms) == 0:
        errmsg = "could not find binding atoms"
        raise RuntimeError(errmsg)

    # Loop over found indexes and see if any of them have unique bonded atoms
    for index in [ix for ix in atg1.atoms.ix]:
        sel_str = (f"(bonded index {index}) and not (index {anchor_ix}) "
                   f"{exclusion_str}")
        atg2 = atomgroup.select_atoms(sel_str)
        if atg2.n_atoms > 0:
            break

    try:
        return index, atg2.atoms[0].ix
    except IndexError:
        errmsg = "could not find third atom"
        raise RuntimeError(errmsg) from None


class FindBindingSite(AnalysisBase):
    """Class to get binding site residues based on the ligand contacts during
    the length of the trajectory

    Parameters
    ----------
    ligand : MDAnaysis.AtomGroup
        Ligand atoms to search for contacts around
    host : MDAnalysis.AtomGroup
        Host atoms (i.e. protein) to search for potential contact residues.
    contact_cutoff : float [5.0]
        Maximum distance in Angstroms to search for contacts.
    contact_precentage : float [20.0]
        Minimum percentage of frames within contact distance to be considered
        as a potential binding site residue. Note, a contact percentage of 0%
        will return all residues in the host.
    residue_selection : string ['backbone']
        Selection string for which atoms to add to the binding site atom group
        from each contact residue.

    Attributes
    ----------
    binding_site : MDAnalysis.AtomGroup
        Binding site atoms as identified from contacts, and sub-selected by the
        `residue_selection` paramter.
    contact_resindices : numpy array
        One dimensional array of the residues found to be in contact of the
        ligand more than contact_percentage of the time.

    Raises
    ------
    UserWarning
        * If fewer than 3 binding site residues have been found.
    """
    def __init__(self, ligand, host, contact_cutoff: float = 5.0,
                 contact_precentage: float = 20.0,
                 residue_selection: str = "backbone", **kwargs):
        super(FindBindingSite, self).__init__(ligand.universe.trajectory,
                                              **kwargs)
        self.ligand = ligand
        self.host = host
        self.contact_cutoff = contact_cutoff
        self.contact_percentage = contact_precentage
        self.residue_selection = residue_selection

    def _prepare(self):
        # accumulator dictionary of length universe residues
        self._host_contacts = {index: 0 for index in
                               self.host.residues.resindices}

    def _single_frame(self):
        pairs = _search_from_capped(self.ligand, self.host,
                                    self.ligand.universe, self.contact_cutoff,
                                    return_distances=False)

        # Get unique atoms in host
        host_atoms = np.unique(pairs[:, 1])
        resixs = self.host.atoms[host_atoms].residues.resindices

        # feed the accumulator
        for entry in resixs:
            self._host_contacts[entry] += 1

    def _conclude(self):
        cutoff_val = (self.n_frames / 100) * self.contact_percentage

        # list of resids with > contact_percentage contacts
        self.contact_resindices = []

        # loop over the residues in the host
        for entry in self._host_contacts:
            if self._host_contacts[entry] >= cutoff_val:
                self.contact_resindices.append(entry)

        # Warn if fewer than 3 residues found
        if len(self.contact_resindices) < 3:
            wmsg = "Fewer than 3 residues identified in binding site"
            warnings.warn(wmsg)

        self.contact_resindices = np.asarray(self.contact_resindices)
        self.binding_site = self.host.residues[
            self.contact_resindices].atoms.select_atoms(self.residue_selection)


class FindHostAtoms(AnalysisBase):
    """Class to get a list of restraint host atoms using a distance search
    around the ligand bonding atoms.

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
        unless `protein_routine` is `False`. ["protein and name CA"]
    num_restraints : int
        Minimum number of boresch restraints to try to generate for a given
        ligand anchor atoms. [3]
    protein_routine : bool
        Option to turn off the C/N atom gathering routine if `p_selection` is
        passed as "protein and name CA". [True]
    search_init_cutoff : float
        Minimum cutoff distance to look for host anchor atoms. [5.0]
    search_max_cutoff : float
        Maximum cutoff distance to look for host anchor atoms. [9.0]

    Notes
    -----
        * The host anchor search will first try to find a minimum of
          ```num_restraints`` anchors within the ``search_init_cutoff``
          distance. Failing to do this the distance will increase by 1 A until
          it exceeds the ``search_max_cutoff`` value. If it fails to find a
          number of host anchors >= ``num_restraints`` it will throw a user
          warning.
    """
    def __init__(self, atomgroup, l_atom, p_selection="protein and name CA",
                 num_restraints=3, protein_routine=True,
                 search_init_cutoff=5, search_max_cutoff=9, **kwargs):
        super(FindHostAtoms, self).__init__(atomgroup.universe.trajectory,
                                            **kwargs)
        self.atomgroup = atomgroup
        self.l_atom = l_atom
        self.p_selection = p_selection
        self.ligand_ag = atomgroup.select_atoms(f'index {l_atom}')
        if len(self.ligand_ag) > 1:
            raise ValueError('Too many ligand atoms passed.')
        self.protein_ag = atomgroup.select_atoms(p_selection)
        self.num_restraints = num_restraints
        self.protein_routine = protein_routine
        self.search_init_cutoff = search_init_cutoff
        self.search_max_cutoff = search_max_cutoff

    def _prepare(self):
        self.host_atoms = []
        self._anchor_pairs = []
        self._anchor_distances = []

    def _single_frame(self):
        pairs, distances = _search_from_capped(self.ligand_ag,
                                               self.protein_ag,
                                               self.atomgroup.universe,
                                               self.search_max_cutoff)
        self._anchor_pairs.append(pairs)
        self._anchor_distances.append(distances)

    def _conclude(self):
        # first we check that we found sufficient numbers of anchors
        found_atoms = 0
        cutoff_distance = self.search_init_cutoff
        while found_atoms < self.num_restraints:
            anchors = []
            for pairs, distances in zip(self._anchor_pairs,
                                        self._anchor_distances):
                for pair, distance in zip(pairs, distances):
                    if distance < cutoff_distance:
                        index = self.protein_ag.atoms[pair[1]].index
                        anchors.append(index)

            anchors = set(anchors)
            found_atoms = len(anchors)

            if found_atoms < self.num_restraints:
                if cutoff_distance < self.search_max_cutoff:
                    cutoff_distance += 1
                    wmsg = (f"Too few anchor atoms found, expanding cutoff "
                            f"to {cutoff_distance}")
                    warnings.warn(wmsg)
                else:
                    wmsg = (f"Too few anchor atoms found, carrying on with "
                            f"{found_atoms} anchors")
                    warnings.warn(wmsg)
                    break

        for entry in anchors:
            if self.protein_routine and (self.p_selection ==
                                         "protein and name CA"):
                ix2, ix3 = _get_bonded_host_cn_atoms(self.atomgroup, entry)
            else:
                ix2, ix3 = _get_bonded_host_atoms(self.atomgroup, entry)

            self.host_atoms.append([entry, ix2, ix3])


def find_host_atoms(atomgroup, l_atom, p_selection="protein and name CA",
                    num_restraints=3, protein_routine=True,
                    search_init_cutoff=5, search_max_cutoff=9):
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
        Minimum number of boresch restraints to try to generate for a given
        ligand anchor atom. [3]
    protein_routine : bool
        Option to turn off the C/N atom gathering routine if `p_selection` is
        passed as "protein and name CA"
    search_init_cutoff : float
        Minimum cutoff distance to look for host anchor atoms. [5.0]
    search_max_cutoff : float
        Maximum cutoff distance to look for host anchor atoms. [9.0]
    """
    print(p_selection)

    host_anchors = _get_host_anchors(atomgroup, l_atom,
                                     anchor_selection=p_selection,
                                     num_atoms=num_restraints,
                                     init_cutoff=search_init_cutoff,
                                     max_cutoff=search_max_cutoff)

    p_atoms = [[i] for i in host_anchors]

    for entry in p_atoms:
        if protein_routine and (p_selection == "protein and name CA"):
            ix2, ix3 = _get_bonded_host_cn_atoms(atomgroup, entry[0])
        else:
            ix2, ix3 = _get_bonded_host_atoms(atomgroup, entry[0])

        entry.extend([ix2, ix3])

    return p_atoms


def find_ligand_atoms(atomgroup, l_selection="resname LIG and not name H*",
                      num_restraints=3, method="RMSF",
                      p_align="protein and name CA"):
    """Gets a list of restraint ligand atoms.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        System for which we want to search for the restraint.
    l_selection : str
        Ligand atom selection ["resname LIG"]
    num_restraints : int
        Number of ligand restraints to attempt to generate [3]
    method : "RMSF"
        Method to use to get ligand anchor atoms, currently only RMSF is
        implemented.
    p_align : str
        Selection string for aligning the trajectory to the host molecule,
        ['protein and name CA'].

    Returns
    -------
    l_atoms : list of 3 int lists
        List containing the potential ligand restraint atoms.
    """

    if method != "RMSF":
        raise NotImplementedError(f"{method} is not implemented yet")
    else:
        return _get_ligand_atoms_rmsf(atomgroup, l_selection, num_restraints,
                                      p_align)


def _get_ligand_atoms_rmsf(atomgroup, l_selection, num_restraints, p_align):
    """Helper function to find potential ligand restraint atoms using RMSF.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
    l_selection : str
        Ligand atom selection.
    num_restraints : int
        Max number of restraints to generate.
    p_align : str
        Selection string to align host.
    """
    # Let's not alter the original trajectory
    copy_u = atomgroup.universe.copy()
    ligand = copy_u.select_atoms(l_selection)
    copy_u.transfer_to_memory()

    # if alignment selection is empty then raise error and let users know
    if len(copy_u.select_atoms(p_align)) == 0:
        errmsg = (f"no atoms matching '{p_align}' found for alignment ",
                  "please review the selection given to p_align")
        raise RuntimeError(errmsg)

    # Align to initial frame first
    prealigner = align.AlignTraj(copy_u, copy_u, select=p_align,
                                 in_memory=True).run()

    # Get the reference frame as the average structure
    ref_coords = copy_u.trajectory.timeseries().mean(axis=1)
    ref = mda.Merge(copy_u.atoms).load_new(ref_coords[:, None, :], order="afc")

    # Align to the average coordiantes
    aligner = align.AlignTraj(copy_u, ref, select=p_align,
                              in_memory=True).run()

    rmsfer = RMSF(ligand).run()

    # Create dictionary of index : rmsf pairs
    rmsf_pairs = {i: x for i, x in zip(ligand.atoms.ix, rmsfer.rmsf)}

    # Get the num_restraints top indices are potential anchor points
    anchors = sorted(rmsf_pairs, key=rmsf_pairs.get)[:num_restraints]

    # loop through anchors and get the lowest RMSF angle set
    l_atoms = []

    for ix in anchors:
        anchor_atg = copy_u.atoms[ix]
        angles = anchor_atg.angles.atomgroup_intersection(ligand, strict=True)
        angle_list = []
        for angle in angles.indices:
            if (angle[0] == ix) or (angle[-1] == ix):
                ixlist = angle.tolist()
                score = (rmsf_pairs[angle[0]] + rmsf_pairs[angle[1]] +
                         rmsf_pairs[angle[2]])
                if ixlist[-1] == ix:
                    ixlist.reverse()
                angle_list.append((ixlist, score))
        angle_list.sort(key=lambda x: x[1])
        l_atoms.append(angle_list[0][0])

    return l_atoms
