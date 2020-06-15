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

    def store(self):
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
        @staticmethod
        def shift_periodic(in_values, periodic=False):
            prob_data = in_values.copy()
            if periodic:
                arange = np.max(prob_data) - np.min(prob_data)
                if arange > 180:
                    for i in range(len(prob_data)):
                        if prob_data[i] < 0:
                            prob_data[i] += 360
            return prob_data

        @staticmethod
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
        titlestring = f"Histogram of {atype} distribution"
        plt.title(titlestring)
        xstring = f"{atype} [{units}]"
        plt.xlabel(xstring)
        plt.ylabel("Number of frames")
        filename = f"{atype}.png"
        plt.savefig(filename, format="png")
        plt.close()


class Bond(VectorData):
    def __init__(self, ix1, ix2, universe):
        """Initialise the Bond object.

        Input
        -----
        ix1 : index of the first atom involved in bond
        ix2 : index of the second atom involve in bond
        universe: MDAnalysis universe instance defining the bond
        """
        # Set values from VectorData
        super(Bond, self).__init__(universe.trajectory.n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2], universe)
        self.atype = "bond"
        self.periodic = False
        self.units = "Å"

    def store(self):
        """Store the current timestep's bond value"""
        index = self.atomgroup.universe.trajectory.frame
        self.values[index] = self.atomgroup.bond.length()


class Angle(VectorData):
    def __init__(self, ix1, ix2, ix3, universe):
        """Initialise the Angle object.

        Input
        -----
        ix1 : index of the first atom involved in angle
        ix2 : index of the second atom involved in angle
        ix3 : index of third atom involved in angle
        universe: MDAnalysis universe instance defining the bond
        """
        # Set value from VectorData
        super(Angle, self).__init__(universe.trajectory.n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2, ix3], universe)
        self.atype = "angle"
        self.periodic = True
        self.units = "degrees"

    def store(self):
        """Store the current timestep's angle value
        """
        index = self.atomgroup.universe.trajectory.frame
        self.values[index] = self.atomgroup.angle.value()


class Dihedral(VectorData):
    def __init__(self, ix1, ix2, ix3, ix4, universe):
        """Initialise the Dihedral object.

        Input
        -----
        ix1 : index of the first atom involved in angle
        ix2 : index of the second atom involved in angle
        ix3 : index of third atom involved in angle
        ix4 : index of the fourth atom involved in angle
        universe: MDAnalysis universe instance defining the bond
        """
        # Set value from VectorData
        super(Dihedral, self).__init__(universe.trajectory.n_frames)
        # We generate the atom group from the zero-based indices
        self.atomgroup = mda.AtomGroup([ix1, ix2, ix3, ix4], universe)a
        self.atype = "dihedral"
        self.periodic = True
        self.units = "degrees"

    def store(self, index):
        """Store the current timestep's dihedral value
        """
        index = self.atomgroup.universe.trajectory.frame
        self.data.values[index] = self.atomgroup.dihedral.value()

