"""
A routine to analyse a trajectory for potential BoreschRestraint parameters.
Dependencies:
    - MDAnalysis
"""

import MDAnalysis as mda
import numpy as np
from scipy import stats
from scipy.stats import circmean, circvar, circstd
from matplotlib import pyplot as plt


class VectorData:
    def __init__(self, n_frames):
        """Initliase the vector data object.

        Input
        -----
        n_frames : int, the number of frames in trajectory
        """
        self.values = np.zeros(n_frames)

    def analyze(self, atype, periodic=True):
        """Analyzes the data held in numpy vector"""
        if not periodic:
            self.mean = self.values.mean()
            self.stdev = self.values.std()
            self.var = self.values.var()
        else:
            p_high = 180
            if atype == "angle":
                p_low = 0
            else:
                p_low = -180
            self.mean = circmean(self.values, low=p_low, high=p_high)
            self.stdev = circstd(self.values, low=p_low, high=p_high)
            self.var = circvar(self.values, low=p_low, high=p_high)

        # Analyze normality
        # Note: this will be broken for distributions that go over a period
        # To be fixed
        self.anderson = stats.anderson(self.values, dist='norm')

    def mean_squared(self):
        """Returns (value-mean)**2 """
        self.ms_values = (self.values - self.mean)**2

    def plot(self, atype, units, periodic=False):
        """Plots the data

        Input
        -----
        atype : str, the data type (e.g. Bond, Dihedral_Angle_1)
        units : units of data (e.g. Å or degrees)

        Notes
        -----
        This module will currently yield odd results for distributions
        that go beyond a period. To be fixed in #2. - IA
        """
        # Set some parameters and prep figure
        pltparams = {'font.size': 12,
                     'axes.titlesize': 16,
                     'axes.titlepad': 15}
        plt.rcParams.update(pltparams)
        plt.figure(1, figsize=[16, 12])

        # Create ProbPlot plot for normality evaluation
        plt.subplot(221)

        # Test for range > 180 degrees and shift for probplots
        if periodic:
            prob_data = self.values.copy()
            arange = np.max(prob_data) - np.min(prob_data)
            if arange > 180:
                for i in range(len(prob_data)):
                    if prob_data[i] < 0:
                        prob_data[i] += 360
        else:
            prob_data = self.values.copy()

        res = stats.probplot(prob_data, dist='norm', plot=plt)

        # Create histogram of quantity
        plt.subplot(222)
        n, bins, patches = plt.hist(self.values, 50, color='c',
                                    edgecolor='k', alpha=0.5)
        # Plot vertical line of mean value
        plt.axvline(self.mean, color='r', linestyle='dashed', linewidth=3)
        # Set plot variables and save to png
        titlestring = f"Histogram of {atype} distribution"
        plt.title(titlestring)
        xstring = f"{atype} [{units}]"
        plt.xlabel(xstring)
        plt.ylabel("Number of frames")
        filename = f"{atype}.png"
        plt.savefig(filename, format="png")
        plt.close()


class BoreschRestraint:
    """A class to store and analyze the bond/angle/diehdral information related
    to a given Boresch restraint."""

    def __init__(self, l_atoms, p_atoms, universe):
        """Init routine for the BoreschRestraint class.
        Input
        -----
        l_atoms : list (size 3) of indices of the ligand atoms
        p_atoms : list (size 3) indices of the protein atoms
        universe : MDAnalysis universe instance

        Note
        ----
        Bond is defined by (l_atoms[0], p_atoms[0])
        Angles are defined by:
            (l_atoms[1], l_atoms[0], p_atoms[0])
            (l_atoms[0], p_atoms[0], p_atoms[1])
        Dihedrals are defined by:
            (l_atoms[2], l_atoms[1], l_atoms[0], p_atoms[0])
            (l_atoms[1], l_atoms[0], p_atoms[0], p_atoms[1])
            (l_atoms[0], p_atoms[0], p_atoms[1], p_atoms[2])
        """
        # Create the Bond, Angle and Dihedral lists
        self.bond = self.Bond(l_atoms[0], p_atoms[0], universe)
        self.angles = []
        self.angles.append(self.Angle(l_atoms[1], l_atoms[0], p_atoms[0],
                                      universe))
        self.angles.append(self.Angle(l_atoms[0], p_atoms[0], p_atoms[1],
                                      universe))
        self.dihedrals = []
        self.dihedrals.append(self.Dihedral(l_atoms[2], l_atoms[1], l_atoms[0],
                                            p_atoms[0], universe))
        self.dihedrals.append(self.Dihedral(l_atoms[1], l_atoms[0], p_atoms[0],
                                            p_atoms[1], universe))
        self.dihedrals.append(self.Dihedral(l_atoms[0], p_atoms[0], p_atoms[1],
                                            p_atoms[2], universe))

    def store_frame(self, index):
        """Function to store data for objects

        Input
        -----
        index : current frame number
        """
        self.bond.store(index)
        self.angles[0].store(index)
        self.angles[1].store(index)
        self.dihedrals[0].store(index)
        self.dihedrals[1].store(index)
        self.dihedrals[2].store(index)

    def analyze(self):
        """Function to analyse and store data"""
        # Deal with bond, angles and dihedrals
        self.bond.data.analyze(atype="bond", periodic=False)
        self.angles[0].data.analyze(atype="angle", periodic=True)
        self.angles[1].data.analyze(atype="angle", periodic=True)
        self.dihedrals[0].data.analyze(atype="dihedral", periodic=True)
        self.dihedrals[1].data.analyze(atype="dihedral", periodic=True)
        self.dihedrals[2].data.analyze(atype="dihedral", periodic=True)
        self.varsum = self.__sum_var()

    def plot(self):
        """Function to plot all the analyzed data"""
        # Deal with bond
        self.bond.data.plot(atype="bond", units="Å", periodic=False)
        # We call them angle1 and angle2
        self.angles[0].data.plot(atype="angle1", units="degrees",
                                 periodic=True)
        self.angles[1].data.plot(atype="angle2", units="degrees",
                                 periodic=True)
        # Similarly dihedral1, dihedral2, dihedral3
        self.dihedrals[0].data.plot(atype="dihedral1", units="degrees",
                                    periodic=True)
        self.dihedrals[1].data.plot(atype="dihedral2", units="degrees",
                                    periodic=True)
        self.dihedrals[2].data.plot(atype="dihedral3", units="degrees",
                                    periodic=True)

    def __sum_var(self):
        """Helper function to add the variances of all varialbles.

        Note: as an initial attempt, we assume all variances to be independent
        (which is known to not be the case). In the future we shall attempt to
        account for the covariance between the bonds/angles.
        """
        combined_variance = (self.bond.data.var +
                             self.angles[0].data.var +
                             self.angles[1].data.var +
                             self.dihedrals[0].data.var +
                             self.dihedrals[1].data.var +
                             self.dihedrals[2].data.var)
        return combined_variance

    def rmsd(self):
        """Helper function to calculate the rmsd of all frames based on the
        mean BoreschRestraint bond/angles.
        """
        self.bond.data.mean_squared()
        self.angles[0].data.mean_squared()
        self.angles[1].data.mean_squared()
        self.dihedrals[0].data.mean_squared()
        self.dihedrals[1].data.mean_squared()
        self.dihedrals[2].data.mean_squared()

        self.rmsd_values = (self.bond.data.ms_values +
                            self.angles[0].data.ms_values +
                            self.angles[1].data.ms_values +
                            self.dihedrals[0].data.ms_values +
                            self.dihedrals[1].data.ms_values +
                            self.dihedrals[2].data.ms_values)

        self.rmsd_values = np.sqrt(self.rmsd_values / 6)

    class Bond:
        def __init__(self, atom1, atom2, universe):
            """Initialise the Bond object.

            Input
            -----
            atom1 : index of the first atom involved in bond
            atom2 : index of the second atom involve in bond
            universe: MDAnalysis universe instance defining the bond
            """
            # Store the atoms in 1-based indices
            self.atom1 = atom1
            self.atom2 = atom2
            # We generate the atom group from the zero-based indices
            self.atomgroup = mda.AtomGroup([atom1-1, atom2-1], universe)
            print(self.atomgroup.atoms[0])
            print(self.atomgroup.atoms[1])
            self.data = VectorData(universe.trajectory.n_frames)

        def store(self, index):
            """Store the current timestep's bond value

            Input
            -----
            index : current frame index
            """
            self.data.values[index] = self.atomgroup.bond.length()

    class Angle:
        def __init__(self, atom1, atom2, atom3, universe):
            """Initialise the Angle object.

            Input
            -----
            atom1 : index of the first atom involved in angle
            atom2 : index of the second atom involved in angle
            atom3 : index of third atom involved in angle
            universe: MDAnalysis universe instance defining the bond
            """
            # Store the atoms in 1-based indices
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            # We generate the atom group from the zero-based indices
            self.atomgroup = mda.AtomGroup([atom1-1, atom2-1, atom3-1],
                                           universe)
            self.data = VectorData(universe.trajectory.n_frames)

        def store(self, index):
            """Store the current timestep's angle value

            Input
            -----
            index : current frame index
            """
            self.data.values[index] = self.atomgroup.angle.value()

    class Dihedral:
        def __init__(self, atom1, atom2, atom3, atom4, universe):
            """Initialise the Dihedral object.

            Input
            -----
            atom1 : index of the first atom involved in angle
            atom2 : index of the second atom involved in angle
            atom3 : index of third atom involved in angle
            atom4 : index of the fourth atom involved in angle
            universe: MDAnalysis universe instance defining the bond
            """
            # Store the atoms in 1-based indices
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.atom4 = atom4
            # We generate the atom group from the zero-based indices
            self.atomgroup = mda.AtomGroup([atom1-1, atom2-1, atom3-1,
                                            atom4-1], universe)
            self.data = VectorData(universe.trajectory.n_frames)

        def store(self, index):
            """Store the current timestep's dihedral value

            Input
            -----
            index : current frame index
            """
            self.data.values[index] = self.atomgroup.dihedral.value()


def __get_alpha_carbons(l_atom, universe, num_atoms=3, max_cutoff=9):
    """Helper function to get the alpha carbons nearest to l_atom.

    This function will first start by finding available atoms within a 5 Å
    cutoff. If < 3 suitable atoms are available, it will expand the cut-off by
    1 Å until it hits a {max_cutoff} Å cut-off range.

    Input
    -----
    l_atom : int, the index of the ligand atom
    universe: mda.Universe object
    """

    found_atoms = 0
    cutoff = 5

    while found_atoms < num_atoms:
        sel_str = f"(protein and name CA) and around {cutoff} bynum {l_atom}"
        alpha_carbons = universe.select_atoms(sel_str)
        found_atoms = alpha_carbons.n_atoms
        if found_atoms < num_atoms:
            if cutoff < max_cutoff:
                cutoff += 1
                warnmsg = f"Too few CAs found, expanding cutoff to {cutoff} Å"
                print(warnmsg)
            else:
                errmsg = "Not enough CA type atoms found"
                raise RuntimeError(errmsg)

    print(f"Number of CA atoms found {alpha_carbons.n_atoms}")

    return alpha_carbons


def __get_CN_atoms(atom, universe):
    """Helper function to get neighbouring atoms of a CA atom

    Input
    -----
    atom : mda.atomgroup of the CA atom of interest
    universe : mda.universe object instance
    """
    resid = atom.resid
    c_atom_str = f"(backbone) and (name C) and (resid {resid})"
    c_atom = universe.select_atoms(c_atom_str)
    n_atom_str = f"(backbone) and (name N) and (resid {resid+1})"
    n_atom = universe.select_atoms(n_atom_str)

    return c_atom.atoms[0].id, n_atom.atoms[0].id


def __RankByVariance(restraints):
    """Helper function which ranks a list of BoreschRestraint objects
    based on their summed variance.

    Note: it's probably not the most efficient, but we are trying to just rank
    out with anything too high in variance.

    Input
    -----
    restraints : list of BoreschRestraint ojbects

    Returns
    -------
    minx_index : the index of the lowest summed variance BoreschRestraints
    """
    var_list = [restraint.varsum for restraint in restraints]

    # Just call list.index(min(list)) that should return the first element it
    # finds and deal with cases where more than one entry has the same value
    min_index = var_list.index(min(var_list))

    return min_index


def __writeBondRestraint(bond, index, force_constant, rfile):
    """Helper function to write out the bond part of the boresch restraint

    Input
    -----
    bond : the bond part of the BoreschRestraint object
    force_constant : the force constant
    rfile : the output restraint file
    """
    rfile.write('[ bonds ]\n')
    rfile.write(';    ai   aj    type  bA         kA       bB    kB\n')
    atom1 = bond.atom1
    atom2 = bond.atom2
    length = bond.data.values[index] / 10.0
    bond_fc = force_constant * 4.184 * 100
    rfile.write(f"{atom1:>6}{atom2:>6}    6    {length:>6.3f}      0.0   "
                f"{length:>6.3f}   {bond_fc:>6.2f}\n")


def __writeAngleRestraint(angles, index, force_constant, rfile):
    """Helper function to write out the angles part of the boresch restraint

    Input
    -----
    angles : the angles part of the BoreschRestraint object
    force_constant : the force constant
    rfile : the output restraint file
    """
    angle_fc = force_constant * 4.184
    rfile.write('[ angles ]\n')
    rfile.write(';   ai    aj    ak    type    thA      '
                'fcA      thB       fcB\n')
    for angle in angles:
        atom1 = angle.atom1
        atom2 = angle.atom2
        atom3 = angle.atom3
        value = angle.data.values[index]
        rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}       1   {value:>6.3f}"
                    f"    0.0   {value:>6.3f}    {angle_fc:>6.2f}\n")


def __writeDihedralRestraint(dihedrals, index, force_constant, rfile):
    """Helper function to write otu the diherals part of the boresch restraint

    Input
    -----
    dihedrals : the dihedrals part of the BoreschRestraint object
    force_constant : the force constant
    rfile : the output restraint file
    """
    angle_fc = force_constant * 4.184
    rfile.write('[ dihedrals ]\n')
    rfile.write(';   ai    aj    ak    al  type    phiA      fcA    '
                'phiB      fcB\n')
    for dihedral in dihedrals:
        atom1 = dihedral.atom1
        atom2 = dihedral.atom2
        atom3 = dihedral.atom3
        atom4 = dihedral.atom4
        value = dihedral.data.values[index]
        rfile.write(f"{atom1:>6}{atom2:>6}{atom3:>6}{atom4:>6}     2    "
                    f"{value:>6.3f}    0.0    {value:>6.3f}   {angle_fc:>6.2f}\n")


def __writeRestraints(restraint_object, index):
    """Helper function to write out the boresch restraint part of a gromacs
    topology file.

    Input
    -----
    restraint_object : BoreschRestraint object
    """

    # For now we set all force constants to 10.0 kcal mol-1
    force_constant = 10.0

    with open('BoreschRestraint.top', 'w') as rfile:
        rfile.write('\n')
        rfile.write('; restraints\n')
        rfile.write('[ intermolecular_interactions ]\n')
        # Write out the bond section
        __writeBondRestraint(restraint_object.bond, index,
                             force_constant, rfile)

        # Write out the angles
        __writeAngleRestraint(restraint_object.angles, index,
                              force_constant, rfile)

        # Write out the dihedrals
        __writeDihedralRestraint(restraint_object.dihedrals, index,
                                 force_constant, rfile)


def __writeClosestFrame(restraint_object, universe):
    """Find the frame that is closer to the mean value.

    Method
    ------
    Get the root mean squared deviation in all metrics and then find the
    frame that is closest to

    Input
    -----
    restraint_object : BoreschRestraint object
    """
    restraint_object.rmsd()

    # Print out the minimum rmsd
    min_rmsd = np.min(restraint_object.rmsd_values)
    print(f"minimum rmsd value: {min_rmsd}")

    # Get the minimum rmsd value index
    index = restraint_object.rmsd_values.argmin()
    print(f"min index: {index}")

    # Output the bond/angle values for user validation
    bond = restraint_object.bond.data.values[index]
    angle1 = restraint_object.angles[0].data.values[index]
    angle2 = restraint_object.angles[1].data.values[index]
    dih1 = restraint_object.dihedrals[0].data.values[index]
    dih2 = restraint_object.dihedrals[1].data.values[index]
    dih3 = restraint_object.dihedrals[2].data.values[index]

    print("Average frame values:")
    print(f"Bond: {bond}, angles: {angle1} {angle2}")
    print(f"Diheds: {dih1} {dih2} {dih3}")

    universe.trajectory[index]
    system = universe.select_atoms('all')
    system.write('ClosestRestraintFrame.gro')

    # We then write out the topology file to have the same
    # Boresch restraints as this frame (this avoids issues
    # with using ligand restraints in gmx)
    __writeRestraints(restraint_object, index)


def FindBoreschRestraints(l_atoms, universe):
    """A function to search through potential Boresch restraints.

    Takes in a set of ligand picked ligand atoms and then finds the best set
    of protein atoms which would constitute a set of BoreschRestraints based
    on the data gathered from a short MD simulation. It will then dump out
    data on the distribution of the bond/angle/dihedral values and a gromacs
    formatted section for use in a .top file.abs

    Input
    -----
    l_atoms : list (size 3) of ligand atom indices. l_atoms[0] will be part of
    the bond with the protein.
    universe : mda.Universe object instance
    """

    # Fetch a group of alpha carbons near the bond atom
    print("Fetching carbon atoms")
    alpha_carbons = __get_alpha_carbons(l_atoms[0], universe, num_atoms=3)

    # Create BoreschRestraint object for each alpha carbon found (add to list)
    restraints = []
    for atom in alpha_carbons.atoms:
        # Get the attached C and N atoms
        c_atom_id, n_atom_id = __get_CN_atoms(atom, universe)

        # Create list and append the atom ids
        p_atoms = []
        p_atoms.append(atom.id)
        p_atoms.append(c_atom_id)
        p_atoms.append(n_atom_id)

        # Create the BoreshRestraint object
        restraints.append(BoreschRestraint(l_atoms, p_atoms, universe))

    # Next we loop through the trajectory and get the bond/angle/dihedrals
    for index, ts in enumerate(universe.trajectory):
        for restraint in restraints:
            restraint.store_frame(index)

    # We then analyze the data
    for restraint in restraints:
        restraint.analyze()

    # Then we rank the restraints based on their stdev
    toprank_index = __RankByVariance(restraints)

    # Create the distribution plots, the gmx restraint entry and closest frame
    restraints[toprank_index].plot()
    __writeClosestFrame(restraints[toprank_index], universe)


if __name__ == "__main__":

    # Deal with argument parsing
    import argparse

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--gro', default='npt_prod1.gro')
        parser.add_argument('--xtc', default='npt_prod1.xtc')
        parser.add_argument('--l1', default=None)
        parser.add_argument('--l2', default=None)
        parser.add_argument('--l3', default=None)
        parser.add_argument('--zeroformat', default=True)
        args = parser.parse_args()
        return args

    print("Parsing arguments")
    args = parse_args()

    print("Creating universe")
    u = mda.Universe(args.gro, args.xtc)

    if args.zeroformat:
        args.l1 = int(args.l1)+1
        args.l2 = int(args.l2)+1
        args.l3 = int(args.l3)+1

    ligand_atoms = [int(args.l1), int(args.l2), int(args.l3)]

    if None in ligand_atoms:
        errmsg = "Missing ligand atoms"
        raise IOError(errmsg)

    FindBoreschRestraints(l_atoms=ligand_atoms, universe=u)
