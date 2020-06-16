"""
Unit and regression test for the MDRestraintsGenerator package.
"""

# Import package, test suite, and other packages as needed
import MDAnalysis as mda
import MDRestraintsGenerator.utils as utils
from .datafiles import T4_TPR, T4_XTC
from numpy.testing import assert_almost_equal
import warnings
import pytest
import sys


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


def test_no_host_anchors(u):
    """Fails to find any host anchors"""

    errmsg = "Too few anchor type atoms found"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.raises(RuntimeError, match=errmsg):
        utils._get_host_anchors(u, l_atom, anchor_selection, num_atoms=1,
                                init_cutoff=1, max_cutoff=3)


def test_cutoff_warning(u):
    """Throws a warning that cutoff is expanding"""

    errmsg = "too few anchor atoms found, expanding cutoff"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.warns(UserWarning, match=errmsg):
        utils._get_host_anchors(u, l_atom, anchor_selection, num_atoms=3,
                                init_cutoff=5, max_cutoff=9)


def test_basic_bonded(u):
    """Basic test for getting a bonded atom"""

    p_atom = 1322
    expected_second_atom = 1320
    expected_third_atom = 1318 

    second_atom, third_atom = utils._get_bonded_atoms(u, p_atom)

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
        utils._get_bonded_atoms(u, p_atom, exclusion_str)
