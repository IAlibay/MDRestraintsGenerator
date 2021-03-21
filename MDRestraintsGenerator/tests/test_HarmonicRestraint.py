"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindHarmonicRestraint
from .datafiles import (T4_TPR, T4_XTC,
                        T4_FB_NDX, T4_H_MDP, T4_FB_OGRO, T4_H_MDP_DEBUG, )
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import filecmp
import pytest


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


def test_basic_regression(tmpdir, u):
    """Standard regressiont test for FindHarmonicRestraint"""
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    finder = search.FindBindingSite(ligand=lig, host=prot)
    finder.run()

    # check we found the right binding site residues
    site_resids = np.array(
        [77,  83,  84,  86,  87,  90,  98,  99, 101, 102, 110, 113, 116, 117,
         120, 128, 132, 152])

    assert_equal(finder.binding_site.residues.resindices, site_resids)

    harmonic_restraint = FindHarmonicRestraint(
            ligand=lig, binding_site=finder.binding_site)

    harmonic_restraint.run()

    # check we got the right standard volume correction
    dG = harmonic_restraint.restraint.standard_state()

    assert_almost_equal(dG, -7.1414, decimal=3)

    # check that the right minimum frame was picked up
    assert harmonic_restraint.restraint.min_frame == 17

    # write out files
    with tmpdir.as_cwd():
        harmonic_restraint.restraint.write()

        closest_gro = mda.Universe('ClosestRestraintFrame.gro')
        closest_gro_ref = mda.Universe(T4_FB_OGRO)

        assert_almost_equal(closest_gro.atoms.positions,
                            closest_gro_ref.atoms.positions)
        assert filecmp.cmp(T4_FB_NDX, 'harmonic_index.ndx')
        assert filecmp.cmp(T4_H_MDP, 'harmonic.mdp')

        # write with debug on
        harmonic_restraint.restraint.write(debug=True)

        assert filecmp.cmp(T4_H_MDP_DEBUG, 'harmonic.mdp')
