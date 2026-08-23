from ase import Atoms
from cluster_generator.selection import get_cluster_indices, get_support_indices


def test_cluster_indices_are_selected_by_element():
    atoms = Atoms(symbols=["Ti", "O", "Co", "Ru"])
    assert get_cluster_indices(atoms, {"Co", "Ru"}) == [2, 3]


def test_support_indices_exclude_cluster_metals():
    atoms = Atoms(symbols=["Ti", "O", "Co", "Ru"])
    assert get_support_indices(atoms, {"Co", "Ru"}) == [0, 1]