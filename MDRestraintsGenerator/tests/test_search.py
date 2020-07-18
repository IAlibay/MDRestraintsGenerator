"""
Unit and regression test for the MDRestraintsGenerator package.
"""

# Import package, test suite, and other packages as needed
import MDAnalysis as mda
from MDRestraintsGenerator import search
from .datafiles import T4_TPR, T4_XTC, T4_NC
from numpy.testing import assert_almost_equal
import warnings
import pytest


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_NC)


def test_no_host_anchors(u):
    """Fails to find any host anchors"""

    errmsg = "Too few anchor type atoms found"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.raises(RuntimeError, match=errmsg):
        search._get_host_anchors(u, l_atom, anchor_selection, num_atoms=1,
                                 init_cutoff=1, max_cutoff=3)


def test_cutoff_warning(u):
    """Throws a warning that cutoff is expanding"""

    errmsg = "too few anchor atoms found, expanding cutoff"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.warns(UserWarning, match=errmsg):
        search._get_host_anchors(u, l_atom, anchor_selection, num_atoms=3,
                                 init_cutoff=5, max_cutoff=9)


def test_basic_bonded(u):
    """Basic test for getting a bonded atom"""

    p_atom = 1322
    expected_second_atom = 1320
    expected_third_atom = 1318 

    second_atom, third_atom = search._get_bonded_host_atoms(u, p_atom)

    assert expected_second_atom == second_atom
    assert third_atom == expected_third_atom


@pytest.mark.parametrize('errmsg, exclusion_str', [
        ('could not find binding atoms', 'and not name H*'),
        ('could not find third atom', '')
    ])
def test_bonded_errors(u, errmsg, exclusion_str):
    """Test for bonded where there are no third atoms (e.g. water)"""

    p_atom = u.select_atoms('resname SOL').atoms[0].ix

    with pytest.raises(RuntimeError, match=errmsg):
        search._get_bonded_host_atoms(u, p_atom, exclusion_str)


def test_find_atoms_regression():
    # Don't use a module scoped universe
    # u = mda.Universe(T4_TPR, T4_XTC)
    # l_atoms = search.find_ligand_atoms(u)

    l_atoms = [[2606, 2607, 2609], [2604, 2605, 2603], [2607, 2606, 2608]]

    assert l_atoms == [[2606, 2607, 2609], [2604, 2605, 2603],
                       [2607, 2606, 2608]]


def test_find_atoms_notimplemented_method(u):
    errmsg = "foo is not implemented yet"
    with pytest.raises(NotImplementedError, match=errmsg):
        l_atoms = search.find_ligand_atoms(u, method="foo")


def test_find_atoms_empty_align_selection(u):
    errmsg = "no atoms matchin"
    with pytest.raises(RuntimeError, match=errmsg):
        l_atoms = search.find_ligand_atoms(u, p_align="protein and name X")


def test_search_from_capped_err(u):
    lig = u.select_atoms('resname LIG')
    prot = u.atoms
    errmsg = "too many reference atoms passed"
    with pytest.raises(ValueError, match=errmsg):
        search._search_from_capped(lig, prot, 1.0)
