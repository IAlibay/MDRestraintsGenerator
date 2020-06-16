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
