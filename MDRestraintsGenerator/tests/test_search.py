"""
Unit and regression test for the MDRestraintsGenerator package.
"""

# Import package, test suite, and other packages as needed
import MDAnalysis as mda
from MDRestraintsGenerator import search
from .datafiles import T4_TPR, T4_XTC, T4_OGRO
import pytest


@pytest.fixture(scope='module')
def u():
    return mda.Universe(T4_TPR, T4_XTC)


def test_no_host_anchors(u):
    """Fails to find any host anchors"""

    errmsg = "Too few anchor type atoms found"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.raises(RuntimeError, match=errmsg):
        search._get_host_anchors(u, l_atom, anchor_selection, num_atoms=1,
                                 init_cutoff=1, max_cutoff=3)


def test_too_many_ligand_atoms_findhostatoms(u):
    """Throws a ValueError if too many ligand atoms have been passed"""

    errmsg = "Too many ligand atoms passed."

    l_atom = "2611 2612"

    with pytest.raises(ValueError, match=errmsg):
        psearch = search.FindHostAtoms(u, l_atom)


def test_no_host_anchors_findhostatoms(u):
    """Throws a warning that too few anchors have been found"""

    errmsg = "Too few anchor atoms found, carrying on with"

    l_atom = 2611

    psearch = search.FindHostAtoms(u, l_atom, search_init_cutoff=1,
                                   search_max_cutoff=3)

    with pytest.warns(UserWarning, match=errmsg):
        psearch.run()


def test_cutoff_warning(u):
    """Throws a warning that cutoff is expanding"""

    errmsg = "too few anchor atoms found, expanding cutoff"

    l_atom = 2611
    anchor_selection = "protein and name CA"

    with pytest.warns(UserWarning, match=errmsg):
        search._get_host_anchors(u, l_atom, anchor_selection, num_atoms=3,
                                 init_cutoff=5, max_cutoff=9)


def test_cutoff_warning_findhostatoms(u):
    """Throws a warning that cutoff is expanding"""

    errmsg = "Too few anchor atoms found, expanding cutoff"

    l_atom = 2611

    psearch = search.FindHostAtoms(u, l_atom)

    with pytest.warns(UserWarning, match=errmsg):
        psearch.run()


def test_findhostatoms_names(u):
    """Checks that the first, second, and third index belong to the right
    atom types"""

    l_atom = 2611

    psearch = search.FindHostAtoms(u, l_atom)

    psearch.run()

    for atoms in psearch.host_atoms:
        assert u.atoms[atoms[0]].name == "CA"
        assert u.select_atoms(f'index {atoms[0]} and bonded index {atoms[1]}')
        assert u.atoms[atoms[1]].name == "C"
        selstring = (f'index {atoms[1]} and (bonded index {atoms[0]}) and '
                     f'(bonded index {atoms[2]})')
        assert u.select_atoms(selstring)
        assert u.atoms[atoms[2]].name == "N"


def test_findhostatoms_manyatoms_err(u):
    l_atom = "2611 2612"
    errmsg = "Too many ligand atoms passed."
    with pytest.raises(ValueError, match=errmsg):
        search.FindHostAtoms(u, l_atom)


@pytest.mark.parametrize('ix1, ix2, ix3', [
    [2560, 2577, 2579], [2581, 2579, 2577]
])
def test_bonded_cn(u, ix1, ix2, ix3):
    """Tests for getting bonded CN atoms"""
    second_atom, third_atom = search._get_bonded_host_cn_atoms(u, ix1)
    assert second_atom == ix2
    assert third_atom == ix3


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


def test_find_atoms_regression(u):
    l_atoms = search.find_ligand_atoms(u)

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


def test_findbindingsite_singleframe():
    u = mda.Universe(T4_OGRO)
    ligand = u.select_atoms('resname LIG')
    host = u.select_atoms('protein')

    bsite = search.FindBindingSite(ligand, host)
    bsite.run()

    assert len(bsite.binding_site) == 72
    assert len(bsite.contact_resindices) == 18


def test_findbindingsite_userwarn():
    u = mda.Universe(T4_OGRO)
    ligand = u.select_atoms('resname LIG')
    host = u.select_atoms('protein')

    bsite = search.FindBindingSite(ligand, host, contact_cutoff=2.15)

    with pytest.warns(UserWarning, match="Fewer than 3"):
        bsite.run()


@pytest.mark.parametrize('perc, atoms, residues', [
    [1, 108, 27], [20, 72, 18]
])
def test_findbindingsite_multiframe_no_min_perc(u, perc, atoms, residues):
    ligand = u.select_atoms('resname LIG')
    host = u.select_atoms('protein')

    bsite = search.FindBindingSite(ligand, host, contact_precentage=perc)
    bsite.run()

    assert len(bsite.binding_site) == atoms
    assert len(bsite.contact_resindices) == residues
