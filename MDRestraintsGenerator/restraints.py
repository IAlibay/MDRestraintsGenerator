"""
restraints.py
A framework for generating restraints for MD simulations

Contains main restraint object classes
"""

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from .datatypes import BoreschRestraint
import numpy as np
import warnings


class FindBoreschRestraint(AnalysisBase):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a Boresch
    restraint from a simulation"""
    def __init__(self, atomgroup, atom_set,
                 force_constant=10.0, **kwargs):
        """Init routine for the BoreschRestraint class.

        Parameters
        ----------
        atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
            Defines the region we want to generate a restraint for.
        atom_set : list of sets
            A list of set of pairs of ligand (l_atoms) and protein (p_atoms)
            atom lists (size 3) with the following order:
            [ (l_atoms, p_atoms).. ]. Note: this is messy
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
        self.atom_set = atom_set
        self.force_constant = force_constant
        self.closest_frame = None

    def _prepare(self):
        """Generates necessary Bond, Angle, Dihedral containers.

        Notes
        -----
        Still in active development
        """
        # Empty container for all the restraints
        self._restraints = []

        for pair in self.atom_set:
            l_atoms = pair[0]
            p_atoms = pair[1]

            self._restraints.append(BoreschRestraint(self.atomgroup,
                                                     l_atoms, p_atoms,
                                                     self.n_frames))

    def _single_frame(self):
        """Loops through trajectory and aggregates necessary data"""
        for restraint in self._restraints:
            restraint.store_frame(self._frame_index)

    def _conclude(self):
        """Analyses and then outputs the best Boresch restraint"""
        for restraint in self._restraints:
            restraint.analyze()

        # Rank restraints based on how much they fluctuate
        # Assign the best restraint
        # Note: to avoid co-linearity we drop anything where the mean angles
        # are < 25 or > 155 degrees (assign it a var of over 9 thousand)
        var_list = []

        for restraint in self._restraints:
            var = restraint.varsum
            for angle in restraint.angles:
                if angle.mean < 25 or angle.mean > 155:
                    var += 9000
            var_list.append(var)

        self.restraint = self._restraints[var_list.index(min(var_list))]

        # Get rmsd & populate min frame/rmsd
        self.restraint.rmsd()
