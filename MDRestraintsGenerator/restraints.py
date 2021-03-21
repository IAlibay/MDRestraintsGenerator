"""
restraints.py
A framework for generating restraints for MD simulations

Contains main restraint object classes
"""

from MDAnalysis.analysis.base import AnalysisBase
from .datatypes import (BoreschRestraint, FlatBottomRestraint,
                        HarmonicRestraint, )


class FindFlatBottomRestraint(AnalysisBase):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a Flat Bottom
    restraint from a simulation.

    Attributes
    ----------
    restraint : :class:`datatypes.FlatBottomRestraint`
        `FlatBottomRestraint` class instance, from which the following can be
        obtained:
          * A plot of the COM distance values via :meth:`plot`
          * Restrained simulation input files via :meth:`write`
          * Standard state correction via :meth:`standard_state`

    """
    def __init__(self, ligand, binding_site, ligand_name="ligand",
                 binding_site_name="binding_site",
                 force_constant=10.0, **kwargs):
        """Init routine for the FindFlatBottomRestraint class.

        Parameters
        ----------
        ligand : MDAnalysis.AtomGroup
            AtomGroup defining the ligand atoms involved in the COM restraint.
        binding_site : MDAnalysis.AtomGroup
            AtomGroup defining the binding site atoms used for the COM
            restraint.
        ligand_name : str
            Name used for the ligand definition in the GMX index file.
            [`ligand`]
        binding_site_name : str
            Name used for the binding site in the GMX index file.
            [`binding_site`]
        force_constant : float
            Force constant of the flat bottom restraint in kcal mol^-1 A^-2
            [10.0]
        """
        super(FindFlatBottomRestraint, self).__init__(
            ligand.universe.trajectory, **kwargs)
        self.ligand_ag = ligand
        self.binding_site_ag = binding_site
        self.ligand_name = ligand_name
        self.binding_site_name = binding_site_name
        self.force_constant = force_constant

    def _prepare(self):
        """Sets up the restraint object"""
        self.restraint = FlatBottomRestraint(
            atomgroup1=self.ligand_ag, atomgroup2=self.binding_site_ag,
            group1_name=self.ligand_name, group2_name=self.binding_site_name,
            n_frames=self.n_frames)

    def _single_frame(self):
        """Aggregates the necessary data for given frame"""
        self.restraint.store_frame(self._frame_index)

    def _conclude(self):
        """Analyses restraint"""
        self.restraint.analyze()


class FindHarmonicRestraint(FindFlatBottomRestraint):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a COM
    harmonic restraint from a simulation.

    Attributes
    ----------
    restraint : :class:`datatypes.HarmonicRestraint`
        `HarmonicRestraint` class instance, from which the following can be
        obtained:
          * A plot of the COM distance values via :meth:`plot`
          * Restrained simulation input files via :meth:`write`
          * Standard state correction via :meth:`standard_state`

    """
    def __init__(self, ligand, binding_site, ligand_name="ligand",
                 binding_site_name="binding_site",
                 force_constant=10.0, **kwargs):
        """Init routine for the FindHarmonicRestraint class.

        Parameters
        ----------
        ligand : MDAnalysis.AtomGroup
            AtomGroup defining the ligand atoms involved in the COM restraint.
        binding_site : MDAnalysis.AtomGroup
            AtomGroup defining the binding site atoms used for the COM
            restraint.
        ligand_name : str
            Name used for the ligand definition in the GMX index file.
            [`ligand`]
        binding_site_name : str
            Name used for the binding site in the GMX index file.
            [`binding_site`]
        force_constant : float
            Force constant of the harmonic restraint in kcal mol^-1 A^-2
            [10.0]
        """
        super(FindHarmonicRestraint, self).__init__(
            ligand, binding_site, ligand_name, binding_site_name,
            force_constant, **kwargs)

    def _prepare(self):
        """Sets up the restraint object"""
        self.restraint = HarmonicRestraint(
            atomgroup1=self.ligand_ag, atomgroup2=self.binding_site_ag,
            group1_name=self.ligand_name, group2_name=self.binding_site_name,
            n_frames=self.n_frames)


class FindBoreschRestraint(AnalysisBase):
    """MDAnalysis.analysis.AnalysisBase derived class to generate a Boresch
    restraint from a simulation

    Attributes
    ----------
    restraint : :class:`datatypes.BoreschRestraint`
        `FlatBottomRestraint` class instance, from which the following can be
        obtained:
          * Plots of the Boresch restraint components via :meth:`plot`
          * Restrained simulation input files via :meth:`write`
          * Standard state correction via :meth:`standard_state`
    """
    def __init__(self, atomgroup, atom_set,
                 force_constant=10.0, **kwargs):
        """Init routine for the FindBoreschRestraint class.

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
