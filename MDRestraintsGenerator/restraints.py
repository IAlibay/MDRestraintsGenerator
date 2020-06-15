"""
restraints.py
A framework for generating restraints for MD simulations

Contains main restraint object classes
"""

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from datatypes import BoreschRestraint
from utils import get_host_atoms, rank_by_variance, get_closest_frame
from io import write_boresch
import numpy as np
import warnings


class FindBoreschRestraint(AnalysisBase):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a Boresch
    restraint from a simulation"""
    def __init__(self, atomgroup, l_atoms=None, p_atoms=None,
                 l_selection="resname LIG",
                 p_selection="protein and name CA", force_constant=10.0,
                 **kwargs):
        """Init routine for the BoreschRestraint class.

        Parameters
        ----------
        atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
            Defines the region we want to generate a restraint for.
        l_atoms : list (size 3) of atomgroups [None]
            Defines ligand atoms.
        p_atoms : list of lists (size 3) of atomgroups [None]
            Defiles the protein atoms to be investigated
        l_selection : str
            Selection string to define which ligand atoms can be chosen for
            consideration as restraint anchor atoms.
            NOTE: Using this is currently unsupported!
        p_selection : str
            Selection string to define which protein atoms can be chosen for
            consideration as restraint anchor atoms. The default
            "protein and name CA" will trigger a special case where alpha
            carbons will be seeked for anchor atoms and the directly bonded C
            and N atoms will be picked as the remainder of the host atoms.
        force_constant : float
            Force constant for the Boresch restraint (kcal/mol) [10.0]

        Notes
        -----
        For the Boresch restraint, the follow are defined:
        Bond is defined by (l_atoms[0], p_atoms[0])
        Angles are defined by:
            (l_atoms[1], l_atoms[0], p_atoms[0])
            (l_atoms[0], p_atoms[0], p_atoms[1])
        Dihedrals are defined by:
            (l_atoms[2], l_atoms[1], l_atoms[0], p_atoms[0])
            (l_atoms[1], l_atoms[0], p_atoms[0], p_atoms[1])
            (l_atoms[0], p_atoms[0], p_atoms[1], p_atoms[2])
        """
        super(FindBoreschRestraint, self).__init__(
                atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        self.l_atoms = l_atoms
        self.p_atoms = p_atoms
        self.l_selection = l_selection
        self.p_selection = p_selection
        self.closest_frame = None

    def _prepare(self):
        """Generates necessary Bond, Angle, Dihedral containers.

        If l_atoms and/or p_atoms are missing, will automatically get them.
        
        Notes
        -----
        Still in active development
        """
        if l_atoms is None:
            errmsg = "undefined ligand input is not currently supported"
            raise ValueError(errmsg)

        if p_atoms is None:
            wmsg = ("no p_atoms selection is passed, will automatically seek "
                    "for host atoms based on {self.p_selection}")
            warnings.warn(wmsg)

            self.p_atoms = get_host_atoms(self.l_atoms[0], self.p_selection,
                                          self.atomgroup, num_atoms=3)

        # Empty container for all the restraints
        self.restraints = []

        for p_atm in p_atoms:
            self.restraints.append(BoreschRestraint(l_atoms, p_atom,
                                                    self.atomgroup,
                                                    self.n_frames))

    def _single_frame(self):
        """Loops through trajectory and aggregates necessary data"""
        for restraint in self.restraints:
            restraint.store_frame(self._frame_index)

    def _conclude(self):
        """Analyses and then outputs the best Boresch restraint"""
        for restraint in self.restraints:
            restraint.analyse()

        # Rank restraints based on how much they fluctuate
        self.toprank_index = rank_by_variance(restraints)

        # Get the closest frame to the average for the topranked restraint
        self.picked_frame = get_cloest_frame(restraints[self.toprank_index])

        # Plot out the statistics of the restraint
        restraints[self.toprank_index].plot(self.picked_frame)

        # Write out the Boresch restraint
        write_boresch(restraints[self.toprank_index], self.picked_frame)

