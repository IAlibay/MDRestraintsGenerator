"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint
from .datafiles import T4_TPR, T4_XTC, T4_OGRO, T4_OTOP
from numpy.testing import assert_almost_equal
import filecmp
import pytest


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


def test_basic_regression(tmpdir, u):
    """Regression test to check we get the same answer"""
    l_atoms = [2611, 2609, 2607]

    p_atoms = search.find_host_atoms(u, l_atoms[0])

    atom_set = [(l_atoms, p) for p in p_atoms]

    find = FindBoreschRestraint(u, atom_set)

    with tmpdir.as_cwd():
        find.run()

        find.restraint.write()
        dG = find.restraint.standard_state()
        
        u_gro = mda.Universe('ClosestRestraintFrame.gro')
        u_gro_ref = mda.Universe(T4_OGRO)

        assert_almost_equal(u_gro.atoms.positions, u_gro_ref.atoms.positions)
        assert filecmp.cmp(T4_OTOP, 'BoreschRestraint.top')
        assert_almost_equal(dG, -6.592, 2)
