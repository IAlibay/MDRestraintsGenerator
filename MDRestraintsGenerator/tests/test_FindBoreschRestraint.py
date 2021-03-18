"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint
from .datafiles import (T4_TPR, T4_XTC, T4_OGRO, T4_OTOP,
                        CB8_TPR, CB8_XTC, CB8_OGRO, CB8_OTOP,)
from numpy.testing import assert_almost_equal, assert_equal
import filecmp
import pytest


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


@pytest.fixture(scope='module')
def u_guest():
    return mda.Universe(CB8_TPR, CB8_XTC)


def test_basic_regression_oldhostsearch(tmpdir, u):
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


def test_hostguest_regression_oldhostsearch(tmpdir, u_guest):
    """Regression test for a Host/Guest system"""
    l_atoms = [155, 157, 165]

    p_atoms = search.find_host_atoms(u_guest, l_atoms[0],
                                     p_selection="resname HST and not name H*")

    atom_set = [(l_atoms, p) for p in p_atoms]

    find = FindBoreschRestraint(u_guest, atom_set)

    with tmpdir.as_cwd():
        find.run()
        find.restraint.write()
        dG = find.restraint.standard_state()

        u_gro = mda.Universe('ClosestRestraintFrame.gro')
        u_gro_ref = mda.Universe(CB8_OGRO)

        assert_almost_equal(u_gro.atoms.positions, u_gro_ref.atoms.positions)
        assert filecmp.cmp(CB8_OTOP, 'BoreschRestraint.top')
        assert_almost_equal(dG, -6.365, 2)


def test_basic_regression(tmpdir, u):
    """Regression test to check we get the same answer"""
    l_atoms = [2611, 2609, 2607]

    psearch = search.FindHostAtoms(u, l_atoms[0])
    psearch.run(stop=1)

    atom_set = [(l_atoms, p) for p in psearch.host_atoms]

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


def test_hostguest_regression(tmpdir, u_guest):
    """Regression test for a Host/Guest system"""
    l_atoms = [155, 157, 165]

    psearch = search.FindHostAtoms(u_guest, l_atoms[0],
                                   p_selection="resname HST and not name H*")
    psearch.run(stop=1)

    atom_set = [(l_atoms, p) for p in psearch.host_atoms]

    find = FindBoreschRestraint(u_guest, atom_set)

    with tmpdir.as_cwd():
        find.run()
        find.restraint.write()
        dG = find.restraint.standard_state()

        u_gro = mda.Universe('ClosestRestraintFrame.gro')
        u_gro_ref = mda.Universe(CB8_OGRO)

        assert_almost_equal(u_gro.atoms.positions, u_gro_ref.atoms.positions)
        assert filecmp.cmp(CB8_OTOP, 'BoreschRestraint.top')
        assert_almost_equal(dG, -6.365, 2)


def test_basic_regression_ligand_search(u):
    """Regression test to check we get the same answer on a ligand search"""

    ligand_atoms = search.find_ligand_atoms(u)

    assert ligand_atoms == [[2606, 2607, 2609], [2604, 2605, 2603],
                            [2607, 2606, 2608]]

    atom_set = []

    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0])
        psearch.run(stop=1)
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

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


def test_basic_regression_ligand_protein_search(u):
    """Regression test to check we get the same answer on a ligand search"""

    ligand_atoms = search.find_ligand_atoms(u)

    assert ligand_atoms == [[2606, 2607, 2609], [2604, 2605, 2603],
                            [2607, 2606, 2608]]

    atom_set = []

    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0])
        psearch.run()
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

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
