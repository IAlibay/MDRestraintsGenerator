"""
Unit and regression test for the MDRestraintsGenerator package.
"""

# Import package, test suite, and other packages as needed
from MDRestraintsGenerator import writers
import pytest
import filecmp
from io import FileIO


@pytest.mark.parametrize('inval, retval', [
    [True, 'yes'], [False, 'no']
])
def test_yes_no(inval, retval):
    assert writers.yes_no(inval) == retval


def test_pull_groups_valueerror(tmpdir):
    with tmpdir.as_cwd():
        with open('test.mdp', 'w') as filed:
            with pytest.raises(ValueError, match="different length"):
                writers._write_pull_groups(filed, ['foo', 'bar'], [12, ])


pull_header = """\
;----------------------------------------------------
; PULL RESTRAINT OPTIONS
;----------------------------------------------------
pull = no
pull-print-com = yes
pull-print-ref-value = yes
pull-print-components = yes
pull-nstxout = 42
pull-nstfout = 42
pull-pbc-ref-prev-step-com = no
pull-xout-average = yes
pull-fout-average = yes
"""


def test_pull_header(tmpdir):
    with tmpdir.as_cwd():
        with open('test.mdp', 'w') as filed:
            writers._write_pull_header(
                filed, pull=False, pull_print_com=True,
                pull_print_ref_value=True, pull_print_components=True,
                pull_nstxout=42, pull_nstfout=42,
                pull_pbc_ref_prev_step_com=False, pull_xout_average=True,
                pull_fout_average=True)

        assert pull_header == open('test.mdp').read()


pull_groups = """\
pull-ngroups = 3
pull-group1-name = foo
pull-group2-name = bar
pull-group3-name = pie
pull-group1-pbcatom = 2
pull-group2-pbcatom = 25
pull-group3-pbcatom = 43
"""


def test_pull_groups(tmpdir):
    with tmpdir.as_cwd():
        with open('test.mdp', 'w') as filed:
            writers._write_pull_groups(filed, group_names=['foo', 'bar', 'pie'],
                                       group_pbc_atoms=[1, 24, 42])

        assert pull_groups == open('test.mdp').read()


pull_coords = """\
pull-ncoords = 2
pull-coord1-type = umbrella
pull-coord1-geometry = distance
pull-coord1-groups = 1 2
pull-coord1-dim = Y Y Y
pull-coord1-start = yes
pull-coord1-init = 0.123
pull-coord1-rate = 0.000
pull-coord1-k = 10.000
pull-coord1-kB = 0.000
pull-coord1-potential-provider = 
pull-coord1-origin = 0.0 0.0 0.0
pull-coord1-vec = 0.0 0.0 0.0
pull-coord2-type = flat-bottom
pull-coord2-geometry = cylinder
pull-coord2-groups = 3 4
pull-coord2-dim = N N Y
pull-coord2-start = no
pull-coord2-init = 0.000
pull-coord2-rate = 1.000
pull-coord2-k = 25.234
pull-coord2-kB = 0.123
pull-coord2-potential-provider = foobar
pull-coord2-origin = 2.56 3.24 0.96
pull-coord2-vec = 0.0 0.0 0.0
"""


def test_pull_coords(tmpdir):
    with tmpdir.as_cwd():
        with open('test.mdp', 'w') as filed:
            writers._write_pull_coords(
                filed,
                pull_coord_types=['umbrella', 'flat-bottom'],
                pull_coord_geometries=['distance', 'cylinder'],
                pull_coord_groups=[(1, 2), (3, 4)],
                pull_coord_dims=['Y Y Y', 'N N Y'],
                pull_coord_starts=[True, False],
                pull_coord_inits=[0.1234, 0.0],
                pull_coord_rates=[0, 1],
                pull_coord_ks=[10, 25.2344],
                pull_coord_kBs=[0.0, 0.1234],
                pull_coord_potential_providers=['', 'foobar'],
                pull_coord_origins=['0.0 0.0 0.0', '2.56 3.24 0.96'],
                pull_coord_vecs=['0.0 0.0 0.0', '0.0 0.0 0.0'])

        assert pull_coords == open('test.mdp').read()
