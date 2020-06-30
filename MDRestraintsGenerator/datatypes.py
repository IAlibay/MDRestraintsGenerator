"""
MDRestraintsGenerator:
A framework for generating restraints for MD simulations

This file contains all the necessary datatypes which are used by the restraint
classes.

Contents:
    - VectorData: template class to store and analyze data along a trajectory.
      - Bond(VectorData)
      - Angle(VectorData)
      - Dihedral(VectorData)

"""

import MDRestraintsGenerator.writers as writers

import MDAnalysis as mda
import numpy as np
from scipy import stats
from scipy.stats import circmean, circvar, circstd
import matplotlib as mpl
from matplotlib import pyplot as plt


class VectorData:
    def __init__(self, n_frames):
        """Initliase the vector data object.

        Input
        -----
        n_frames : int, the number of frames in trajectory
        """
        self.values = np.zeros(n_frames)

    def store(self, index):
        """Store the current timestep's value"""
        raise NotImplementedError("Only implemented in child classes")

    def analyze(self):
        """Analyzes the data held in numpy vector"""
        if not self.periodic:
            self.mean = self.values.mean()
            self.stdev = self.values.std()
            self.var = self.values.var()
        else:
            p_high = 180
            if self.atype == "angle":
                p_low = 0
            else:
                p_low = -180
            self.mean = circmean(self.values, low=p_low, high=p_high)
            self.stdev = circstd(self.values, low=p_low, high=p_high)
            self.var = circvar(self.values, low=p_low, high=p_high)

    def mean_squared(self):
        """Returns (value-mean)**2 """
        self.ms_values = (self.values - self.mean)**2

    def plot(self, picked_frame=None):
        """Plots the data

        Input
        -----
        picked_frame : int, index of picked frame [None]

        Notes
        -----
        This module will currently yield odd results for distributions
        that go beyond a period. To be fixed in #2. - IA

        This is definitely not meant to stay here, long term this will
        get moved to a restraint specific method
        """
        def shift_periodic(in_values, periodic=False):
            prob_data = in_values.copy()
            if periodic:
                arange = np.max(prob_data) - np.min(prob_data)
                if arange > 180:
                    for i in range(len(prob_data)):
                        if prob_data[i] < 0:
                            prob_data[i] += 360
            return prob_data

        def get_hist_binrange(atype, in_values):
            if atype == "bond":
                bin_width = 0.2
                return np.arange(in_values.min(), in_values.max() + bin_width,
                                 bin_width)
            if atype == "angle":
                bin_width = 2.0
                return np.arange(0, 180 + bin_width, bin_width)

            if atype == "dihedral":
                bin_width = 4
                return np.arange(-180, 180 + bin_width, bin_width)


        try:
            self.mean
        except AttributeError:
            raise AttributeError("vector object has not been analyzed yet")



        # Set some parameters and prep figure
        pltparams = {'font.size': 12,
                     'axes.titlesize': 16,
                     'axes.titlepad': 15}
        plt.rcParams.update(pltparams)
        plt.figure(1, figsize=[16, 12])

        # Create ProbPlot plot for normality evaluation
        plt.subplot(221)
        prob_vals = shift_periodic(self.values, self.periodic)
        res = stats.probplot(prob_vals, dist='norm', plot=plt)

        # Create histogram of quantity
        plt.subplot(222)
        bin_range = get_hist_binrange(self.atype, self.values)
        n, bins, patches = plt.hist(self.values, color='c', bins=bin_range,
                                    edgecolor='k', alpha=0.5, label='observed')

        # Plot vertical line of mean value
        plt.axvline(self.mean, color='r', linestyle='dashed', linewidth=3,
                    label='mean')

        # Plot the picked frame if passed
        if picked_frame is not None:
            plt.axvline(self.values[picked_frame], color='k',
                        linestyle='dashed', linewidth=3,
                        label='picked')

        # Set plot variables and save to png
        titlestring = f"Histogram of {self.atype} distribution"
        plt.title(titlestring)
        xstring = f"{self.atype} [{self.units}]"
        plt.xlabel(xstring)
        plt.ylabel("Number of frames")
        filename = f"{self.filename}.png"
        plt.legend(loc="best")
        plt.savefig(filename, format="png")
        plt.close()


class Bond(VectorData):
    def __init__(self, ix1, ix2, atomgroup, n_frames, suffix_index=1):
        """Initialise the Bond object.

        Input
        -----
        ix1 : int
            index of the first atom involved in bond
        ix2 : int
            index of the second atom involve in bond
        atomgroup: MDAnalysis.AtomGroup or MDAnalysis.Universe
            MDAnalysis object containing the bond
        n_frames : int
            number of frames restraint object will store
        suffix_index : int
            number to add as a filename index when plotting [1]
        """
        # Set values from VectorData
        super(Bond, self).__init__(n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2], atomgroup.universe)
        self.atype = "bond"
        self.periodic = False
        self.units = "Å"
        self.filename = f"{self.atype}_{suffix_index}"

    def store(self, index):
        """Store the current timestep's bond value"""
        self.values[index] = self.atomgroup.bond.length()


class Angle(VectorData):
    def __init__(self, ix1, ix2, ix3, atomgroup, n_frames, suffix_index=1):
        """Initialise the Angle object.

        Input
        -----
        ix1 : int
            index of the first atom involved in angle
        ix2 : int
            index of the second atom involved in angle
        ix3 : int
            index of third atom involved in angle
        atomgroup : MDAnalysis.AtomGroup or MDAnalysis.Universe
            MDAnalysis object containing the angle
        n_frames : int
            number of frames restraint object will store
        suffix_index : int
            number to add as a filename index when plotting [1]
        """
        # Set value from VectorData
        super(Angle, self).__init__(n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2, ix3], atomgroup.universe)
        self.atype = "angle"
        self.periodic = True
        self.units = "degrees"
        self.filename = f"{self.atype}_{suffix_index}"

    def store(self, index):
        """Store the current timestep's angle value
        """
        self.values[index] = self.atomgroup.angle.value()


class Dihedral(VectorData):
    def __init__(self, ix1, ix2, ix3, ix4, atomgroup, n_frames,
                 suffix_index=1):
        """Initialise the Dihedral object.

        Input
        -----
        ix1 : int
            index of the first atom involved in dihedral
        ix2 : int
            index of the second atom involved in dihedral
        ix3 : int
            index of third atom involved in dihedral
        ix4 : int
            index of the fourth atom involved in dihedral
        atomgroup : MDAnalysis.AtomGroup or MDAnalysis.Universe
            MDAnalysis object containing the dihedral
        n_frames : int
            number of frames restraint object will store
        suffix_index : int
            number to add as a filename index when plotting [1]
        """
        # Set value from VectorData
        super(Dihedral, self).__init__(n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2, ix3, ix4],
                                       atomgroup.universe)
        self.atype = "dihedral"
        self.periodic = True
        self.units = "degrees"
        self.filename = f"{self.atype}_{suffix_index}"

    def store(self, index):
        """Store the current timestep's dihedral value
        """
        self.values[index] = self.atomgroup.dihedral.value()


class BoreschRestraint:
    """A class to store and analyze the bond/angle/diehdral information related
    to a given Boresch restraint."""
    def __init__(self, atomgroup, l_atoms, p_atoms, n_frames=None):
        """Init routine for the BoreschRestraint class.

        Parameters
        ----------
        atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
            System for which we want to search for the restraint.
        l_atoms : list (size 3) of ints
            Indices of the ligand atoms.
        p_atoms : list (size 3) of int
            Indices of the protein atoms.
        n_frames : int [`None`]
            Number of frames to store data for.
        
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
        # Store atomgroup for later use
        # call atoms in case its a Universe
        self.atomgroup = atomgroup.atoms

        # Either default to all frames or set subset of frames
        if n_frames is None:
            self.n_frames = atomgroup.universe.trajectory.n_frames
        else:
            self.n_frames = n_frames

        # Create the bond, angle and dihedral objects
        self.bond = Bond(l_atoms[0], p_atoms[0], atomgroup, self.n_frames)
        self.angles = []
        self.angles.append(Angle(l_atoms[1], l_atoms[0], p_atoms[0],
                                 atomgroup, self.n_frames, suffix_index=1))
        self.angles.append(Angle(l_atoms[0], p_atoms[0], p_atoms[1],
                                 atomgroup, self.n_frames, suffix_index=2))
        self.dihedrals = []
        self.dihedrals.append(Dihedral(l_atoms[2], l_atoms[1], l_atoms[0],
                                       p_atoms[0], atomgroup, self.n_frames,
                                       suffix_index=1))
        self.dihedrals.append(Dihedral(l_atoms[1], l_atoms[0], p_atoms[0],
                                       p_atoms[1], atomgroup, self.n_frames,
                                       suffix_index=2))
        self.dihedrals.append(Dihedral(l_atoms[0], p_atoms[0], p_atoms[1],
                                       p_atoms[2], atomgroup, self.n_frames,
                                       suffix_index=3))

    def store_frame(self, index):
        """Function to store data for objects

        Paramters
        ---------
        index : current frame number
        """
        self.bond.store(index)

        for angle in self.angles:
            angle.store(index)

        for dihedral in self.dihedrals:
            dihedral.store(index)

    def analyze(self):
        """Function to analyze boresch restraint object data"""
        self.bond.analyze()

        for angle in self.angles:
            angle.analyze()

        for dihedral in self.dihedrals:
            dihedral.analyze()

        self.varsum = self._sum_var()

    def _sum_var(self):
        """Helper function to add the variances of all varialbles.

        Note: as an initial attempt, we assume all variances to be independent
        (which is known to not be the case). In the future we shall attempt to
        account for the covariance between the bonds/angles.
        """
        combined_variance = (self.bond.var +
                             self.angles[0].var +
                             self.angles[1].var +
                             self.dihedrals[0].var +
                             self.dihedrals[1].var +
                             self.dihedrals[2].var)
        return combined_variance

    def rmsd(self):
        """Helper function to calculate the rmsd of all frames based on the
        mean BoreschRestraint bond/angles.

        This generates min_frame and min_rmsd for later use.
        """
        self.bond.mean_squared()
        self.angles[0].mean_squared()
        self.angles[1].mean_squared()
        self.dihedrals[0].mean_squared()
        self.dihedrals[1].mean_squared()
        self.dihedrals[2].mean_squared()

        self.rmsd_values = (self.bond.ms_values +
                            self.angles[0].ms_values +
                            self.angles[1].ms_values +
                            self.dihedrals[0].ms_values +
                            self.dihedrals[1].ms_values +
                            self.dihedrals[2].ms_values)

        self.rmsd_values = np.sqrt(self.rmsd_values / 6)

        self.min_frame = self.rmsd_values.argmin()
        self.min_rmsd = np.min(self.rmsd_values)

    def plot(self, frame=None, path=None):
        """Plots all the analyzed data

        Input
        -----
        frame : int
            index of chosen frame to plot.
        path : str
            path to location where files should be written.
        """

        if frame is None:
            try:
                frame = self.min_frame
            except AttributeError:
                pass

        self.bond.plot(picked_frame=frame)
        self.angles[0].plot(picked_frame=frame)
        self.angles[1].plot(picked_frame=frame)
        self.dihedrals[0].plot(picked_frame=frame)
        self.dihedrals[1].plot(picked_frame=frame)
        self.dihedrals[2].plot(picked_frame=frame)

    def write(self, frame=None, path=None, force_constant=10.0, outtype="GMX"):
        """Writes out boresch restraint

        Input
        -----
        frame : int
            index of frame to write out, will default to frame closes to mean.
        path : str
            path to location where files should be written.
        force_constant : float
            strength of the Boresch restraint [10.0 kcal/mol]
        outtype : str
            type of restraint to write, for now only "GMX" is accepted.
        """

        if frame is None:
            try:
                frame = self.min_frame
            except AttributeError:
                raise RuntimeError("no frame defined for writing")

        if outtype is not "GMX":
            raise RuntimeError(f"{outtype} not implemented yet")

        # Final check for co-linearity
        for angle in self.angles:
            if (angle.values[frame] < 25) or (angle.values[frame] > 155):
                errmsg = (f"picked frame contains angle near colinearity ",
                          f"value: {angle.value[frame]}\n",
                          f"This is a bad idea, choose another set of ",
                          f"restraint atoms.")
                raise RuntimeError(errmsg)

        self._write_gmx(index=frame, force_constant=force_constant)

    def _write_gmx(self, index, force_constant):
        # seek corre
        self.atomgroup.universe.trajectory[index]
        self.atomgroup.write('ClosestRestraintFrame.gro')

        # write out top file
        with open('BoreschRestraint.top', 'w') as rfile:
            rfile.write('\n')
            rfile.write('; restraints\n')
            rfile.write('[ intermolecular_interactions ]\n')

            # bond
            writers._write_bond_header(rfile)
            writers._write_bond(self.bond, index, force_constant, rfile)

            # angle
            writers._write_angle_header(rfile)
            for angle in self.angles:
                writers._write_angle(angle, index, force_constant, rfile)

            # dihedral
            writers._write_dihedral_header(rfile)
            for dihedral in self.dihedrals:
                writers._write_dihedral(dihedral, index, force_constant, rfile)

    def standard_state(self, frame=None, force_constant=10.0,
                       temperature=298.15, calc_type="analytical"):
        """Reports the dG_off standard state correction energy for the
        given frame.

        Input
        -----
        frame : int
            index of frame to write out, will default to frame closest to mean.
        force_constant : float
            strength of the Boresch restraint [10.0 kcal/mol]
        temperature : float
            system temperature in Kelvins [298.15]
        calc_type : str
            calculation type, currently only the analytical correction is
            supported but in the future we'd like to add the numerical
            approaches used by Yank.
        """
        if frame is None:
            try:
                frame = self.min_frame
            except AttributeError:
                raise RuntimeError("no frame defined to get energy for")


        if calc_type is not "analytical":
            raise NotImplementedError(f"{calc_type} is not implemented")
        else:
            return self._analytical_energy(frame, force_constant, temperature)

    def _analytical_energy(self, frame, force_constant, temperature):
        """Get the dG_off standard state correction via the Boresch
        analytical method.

        Acknowledgement
        ---------------
        Based on an original script by Matteo Aldeghi.

        Notes
        -----
        This is known to be inaccurate (see Yank).
        """
        Gas_K = (8.314472*0.001) / 4.184 # Gas constant kcal/mol/K
        StandardV = 1.66 # standard volume in nm^3

        rAa = self.bond.values[frame] / 10 # convert to nm
        thA = np.radians(self.angles[0].values[frame])
        thB = np.radians(self.angles[1].values[frame])

        frc_bond = force_constant * 100
        frc_angle = force_constant
        frc_dihe = force_constant

        numerator = 8.0 * (np.pi**2) * StandardV
        numerator *= ((frc_bond * (frc_angle ** 2) * (frc_dihe ** 3)) ** 0.5)
        denominator = (rAa **2) * np.sin(thA) * np.sin(thB)
        denominator *= ((2 * np.pi * Gas_K * temperature) ** 3)

        dG = - Gas_K * temperature * np.log(numerator/denominator)

        return dG

