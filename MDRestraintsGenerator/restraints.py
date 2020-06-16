"""
restraints.py
A framework for generating restraints for MD simulations

Contains main restraint object classes
"""

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from datatypes import BoreschRestraint
from utils import get_host_atoms
from io import write_boresch
import numpy as np
import warnings


class FindBoreschRestraint(AnalysisBase):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a Boresch
    restraint from a simulation"""
    def __init__(self, atomgroup, l_atoms=None, p_atoms=None,
                 l_selection="resname LIG",
                 p_selection="protein and name CA", force_constant=10.0,
                 protein_routine=True, search_init_cutoff=5.0,
                 search_max_cutoff=9.0,
                 **kwargs):
        """Init routine for the BoreschRestraint class.

        Parameters
        ----------
        atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
            Defines the region we want to generate a restraint for.
        l_atoms : list (size 3) of atomgroups [None]
            Defines ligand atoms.
        p_atoms : list of lists (size 3) of atomgroups [None]
            Defines the protein atoms to be investigated
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
        protein_routine : bool
            Option to turn off the C/N atom gathering routine if `p_selection`
            is passed as "protein and name CA".
        search_init_cutoff : float
            Minimum cutoff distance to look for host anchor atoms. Used if
            p_atoms is `None`. [5.0]
        search_max_cutoff : float
            Maximum cutoff distance to look for host anchor atoms. Used if
            p_atoms is `None`. [9.0]

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
        self.force_constant = force_constant
        self.protein_routine = protein_routine
        self.search_init_cutoff = search_init_cutoff
        self.search_max_cutoff = search_max_cutoff
        self.closest_frame = None

    def _prepare(self):
        """Generates necessary Bond, Angle, Dihedral containers.

        If l_atoms and/or p_atoms are missing, will automatically get them.
        
        Notes
        -----
        Still in active development
        """
        if self.l_atoms is None:
            errmsg = "undefined ligand input is not currently supported"
            raise ValueError(errmsg)

        if self.p_atoms is None:
            wmsg = ("no p_atoms selection is passed, will automatically seek "
                    "for host atoms based on {self.p_selection}")
            warnings.warn(wmsg)

            self.p_atoms = get_host_atoms(self.atomgroup, self.l_atoms[0],
                    self.p_selection, num_restraints=3,
                    protein_routine=self.protein_routine,
                    search_init_cutoff=self.search_init_cutoff,
                    search_max_cutoff=self.search_max_cutoff)

        # Empty container for all the restraints
        self.restraints = []

        for p_atm in self.p_atoms:
            self.restraints.append(BoreschRestraint(self.atomgroup,
                                                    self.l_atoms, p_atm,
                                                    self.n_frames))

    def _single_frame(self):
        """Loops through trajectory and aggregates necessary data"""
        for restraint in self.restraints:
            restraint.store_frame(self._frame_index)

    def _conclude(self):
        """Analyses and then outputs the best Boresch restraint"""
        for restraint in self.restraints:
            restraint.analyze()

        # Rank restraints based on how much they fluctuate
        var_list = [restraint.varsum for restrain in restraints]
        self.best_restraint_index = var_lis.index(min(varlist))

        # Get the closest frame to the average for the topranked restraint
        # Get rmsd of all frames from mean
        self.min_rmsd, self.min_frame = _get_min_frame(restraints,
                self.best_restraint_index)

        # Plot out the statistics of the restraint
        restraints[self.toprank_index].plot(self.min_frame)

        # Write out the Boresch restraint
        write_boresch(restraints[self.toprank_index], self.picked_frame)

    @staticmethod
    def _get_min_frame(restraints, index):
        restraints[index].rmsd()
        min_rmsd = np.min(restraints[index].rmsd_values)
        frame = restraints[index].rmsd_values.argmin()

        return min_rmsd, frame
