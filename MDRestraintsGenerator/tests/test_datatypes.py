"""
Unit and regression test for the MDRestraintsGenerator package.
"""

import MDAnalysis as mda
import MDRestraintsGenerator.datatypes as dtypes
from .datafiles import T4_TPR, T4_XTC, T4_OGRO
import pytest
import numpy as np
import os
from pathlib import Path
from numpy.testing import assert_almost_equal


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


@pytest.fixture(scope='module')
def u_gro():
    return mda.Universe(T4_OGRO)


def test_vector_store():
    n_frames = 100
    vector = dtypes.VectorData(n_frames)

    errmsg = "Only implemented in child classes"

    with pytest.raises(NotImplementedError, match=errmsg):
        vector.store(0)


@pytest.mark.parametrize('class_name', ['COGDistance', 'COMDistance'])
def test_cog_com_distance_oneframe_basic(u_gro, class_name):
    ag1 = u_gro.select_atoms('index 0')
    ag2 = u_gro.select_atoms('index 1')
    distance_class = getattr(dtypes, class_name)
    cog_or_com = distance_class(ag1, ag2, n_frames=1)
    cog_or_com.store(0)

    assert_almost_equal(cog_or_com.values, np.array([1.0082664]))


@pytest.mark.parametrize('class_name', ['COGDistance', 'COMDistance'])
def test_cog_com_distance_multiframe_basic(u, class_name):
    ag1 = u.select_atoms('index 0')
    ag2 = u.select_atoms('index 100')
    distance_class = getattr(dtypes, class_name)
    obj = distance_class(ag1, ag2, n_frames=u.trajectory.n_frames)

    for i, ts in enumerate(u.trajectory):
        obj.store(i)

    diff_array = np.array([
        10.873642, 11.237176, 11.031385, 10.970623, 11.89015, 10.771446,
        10.005971, 9.920996, 11.99379, 10.9161215, 11.037367, 9.546233,
        11.420322, 5.6185527, 10.2656975, 9.249049, 10.485258, 8.989676,
        11.3187685, 12.229224, 11.093904], dtype=np.float32)

    assert_almost_equal(obj.values, diff_array, decimal=5)


def test_BaseCOMDistanceRestraint_analyze_error(u):
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    obj = dtypes.BaseCOMDistanceRestraint(lig, prot)

    with pytest.raises(NotImplementedError, match="child"):
        obj.analyze()


@pytest.mark.parametrize('frames', [None, 10])
def test_FlatBottomRestraint_frame(u, frames):
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    if frames is None:
        frames = u.trajectory.n_frames
        obj = dtypes.FlatBottomRestraint(lig, prot)
    else:
        obj = dtypes.FlatBottomRestraint(lig, prot, n_frames=frames)

    assert obj.n_frames == frames
    assert obj.com.values.size == frames


@pytest.mark.parametrize('frame', [None, 10])
def test_FlatBottom_plotting(tmpdir, u, frame):
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.analyze()

    with tmpdir.as_cwd():
        obj.plot(frame=frame)
        assert os.path.isfile('COM-distance_1.png')


def test_FlatBottom_plotting_path(tmpdir, u):
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.analyze()

    with tmpdir.as_cwd():
        obj.plot(frame=None, path='testdir')
        assert os.path.isfile('./testdir/COM-distance_1.png')


def test_FlatBottom_plotting_notanalysed(tmpdir, u):
    lig = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)

    with pytest.raises(AttributeError, match="not been analyzed yet"):
        obj.plot()


def test_FlatBottom_write_path(tmpdir, u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.store_frame(0)
    obj.analyze()

    with tmpdir.as_cwd():
        Path('./testdir').mkdir()
        obj.write(frame=0, path='./testdir')

        u2_gro = mda.Universe('./testdir/ClosestRestraintFrame.gro')
        assert_almost_equal(u_gro.atoms.positions, u2_gro.atoms.positions)
        os.path.isfile('./testdir/flatbottom_index.ndx')


def test_FlatBottom_write_noframe(u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)

    with pytest.raises(RuntimeError, match="no frame defined for writing"):
        obj.write()


def test_FlatBottom_write_no_walldistance(u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.store_frame(0)

    with pytest.raises(RuntimeError, match="The `wall_distance` attribute"):
        obj.write(frame=0)


def test_FlatBottom_write_noimplemented(u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.store_frame(0)
    obj.analyze()

    with pytest.raises(RuntimeError, match="not implemented yet"):
        obj.write(outtype="AMBER")


def test_nearest_com_atom(u_gro):
    ag = u_gro.select_atoms('resid 0')

    ix = dtypes.FlatBottomRestraint._get_nearest_com_atom(ag)

    assert ix == 9


def test_standard_state(u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)
    obj.store_frame(0)

    dG = obj.standard_state(wall_distance=3.3)
    assert_almost_equal(dG, -1.42, decimal=2)


def test_standard_state_no_walldistance_error(u_gro):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    obj = dtypes.FlatBottomRestraint(lig, prot)

    with pytest.raises(AttributeError, match="`wall_distance` was not"):
        obj.standard_state()


@pytest.mark.parametrize('class_name',
                         ['FlatBottomRestraint', 'HarmonicRestraint'])
def test_standard_state_not_implemented_error(u_gro, class_name):
    lig = u_gro.select_atoms('resname LIG')
    prot = u_gro.select_atoms('protein')

    tested_class = getattr(dtypes, class_name)
    obj = tested_class(lig, prot)

    with pytest.raises(NotImplementedError, match="not implemented"):
        obj.standard_state(3.3, calc_type="numerical")


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


def test_Boresch_plotting_path(tmpdir, u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    boresch.analyze()

    with tmpdir.as_cwd():
        boresch.plot(frame=None, path='testdir')
        for name in ['./testdir/bond_1.png', './testdir/angle_1.png',
                     './testdir/angle_2.png', './testdir/dihedral_1.png',
                     './testdir/dihedral_2.png', './testdir/dihedral_3.png']:
            assert os.path.isfile(name)


def test_Boresch_plotting_notanalysed(tmpdir, u):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u, l_atoms, p_atoms)

    errmsg = "vector object has not been analyzed yet"

    with pytest.raises(AttributeError, match=errmsg):
        boresch.plot()


def test_Boresch_write_path(tmpdir, u_gro):
    l_atoms = [0, 1, 2]
    p_atoms = [4, 5, 6]

    boresch = dtypes.BoreschRestraint(u_gro, l_atoms, p_atoms)

    boresch.store_frame(0)

    with tmpdir.as_cwd():
        Path('./testdir').mkdir()
        boresch.write(frame=0, path='./testdir')

        u2_gro = mda.Universe('./testdir/ClosestRestraintFrame.gro')

        assert_almost_equal(u_gro.atoms.positions, u2_gro.atoms.positions)


def test_Boresch_write_colinear(u_gro):
    l_atoms = [2606, 2610, 2611]
    p_atoms = [1218, 1216, 1214]

    boresch = dtypes.BoreschRestraint(u_gro, l_atoms, p_atoms)

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
