"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint
from .datafiles import T4_TPR, T4_XTC, T4_OGRO, T4_OTOP
from MDAnalysisTests.datafiles import GRO, XTC
from numpy.testing import assert_almost_equal, assert_equal
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


def test_aligntraj(tmpdir):
    """AlignTraj is failing, so let's test it here"""
    copy_u = mda.Universe(GRO, XTC)
    with tmpdir.as_cwd():
        prealigner = align.AlignTraj(copy_u, copy_u, select="protein and name CA",
                                     in_memory=True)
    prealigner.run()

    assert 1 == 1


def test_basic_regression_ligand_search(u):
    """Regression test to check we get the same answer on a ligand search"""

    ligand_atoms = search.find_ligand_atoms(u)

    atom_set = []

    for l_atoms in ligand_atoms:
        p_atoms = search.find_host_atoms(u, l_atoms[0])
        atom_set.extend([(l_atoms, p) for p in p_atoms])

    boresch = FindBoreschRestraint(u, atom_set)

    boresch.run()

    assert_equal(boresch.restraint.bond.atomgroup.atoms.ix,
                 [2606, 1563])
    assert_equal(boresch.restraint.angles[0].atomgroup.atoms.ix,
                 [2607, 2606, 1563])
    assert_equal(boresch.restraint.angles[1].atomgroup.atoms.ix,
                 [2606, 1563, 1569])
    assert_equal(boresch.restraint.dihedrals[0].atomgroup.atoms.ix,
                 [2609, 2607, 2606, 1563])
    assert_equal(boresch.restraint.dihedrals[1].atomgroup.atoms.ix,
                 [2607, 2606, 1563, 1569])
    assert_equal(boresch.restraint.dihedrals[2].atomgroup.atoms.ix,
                 [2606, 1563, 1569, 1571])
