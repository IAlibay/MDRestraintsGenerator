"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
import MDRestraintsGenerator.datatypes as dtypes
from .datafiles import T4_TPR, T4_XTC, T4_OGRO
import pytest
import os


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


def test_vector_store():
    n_frames = 100
    vector = dtypes.VectorData(n_frames)

    errmsg = "Only implemented in child classes"

    with pytest.raises(NotImplementedError, match=errmsg):
        vector.store(0)


@pytest.mark.parametrize('frames', [None, 10])
def test_BoreschRestraint_frame(u, frames):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    if frames is None:
        frames = u.trajectory.n_frames
        boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)
    else:
        boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms, frames)

    assert boresch.n_frames == frames

    assert boresch.bond.values.size == frames

    for angle in boresch.angles:
        assert angle.values.size == frames

    for dihed in boresch.dihedrals:
        assert dihed.values.size == frames


@pytest.mark.parametrize('frame', [None, 10])
def test_Boresch_plotting(tmpdir, u, frame):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    boresch.analyze()

    with tmpdir.as_cwd():
        boresch.plot(frame=frame)
        for name in ['bond_1.png', 'angle_1.png', 'angle_2.png',
                     'dihedral_1.png', 'dihedral_2.png', 'dihedral_3.png']:
            assert os.path.isfile(name)


def test_Boresch_plotting_notanalysed(tmpdir, u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "vector object has not been analyzed yet"

    with pytest.raises(AttributeError, match=errmsg):
        boresch.plot()


def test_Boresch_write_colinear():
    u = mda.Universe(T4_OGRO)

    l_atoms = [2606, 2610, 2611]
    p_atoms = [1218, 1216, 1214]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    boresch.store_frame(0)

    with pytest.raises(RuntimeError, match="colinearity"):
        boresch.write(frame=0)


def test_Boresch_write_noframe(u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "no frame defined for writing"

    with pytest.raises(RuntimeError, match=errmsg):
        boresch.write()


def test_Boresch_write_wrongtype(u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "not implemented yet"

    with pytest.raises(RuntimeError, match=errmsg):
        boresch.write(frame=0, outtype="AMBER")


def test_standard_state_noframe_error(u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "no frame defined to get energy for"

    with pytest.raises(RuntimeError, match=errmsg):
        boresch.standard_state()


def test_standard_state_wrongmethod_error(u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "foo is not implemented"

    with pytest.raises(NotImplementedError, match=errmsg):
        boresch.standard_state(frame=0, calc_type='foo')
