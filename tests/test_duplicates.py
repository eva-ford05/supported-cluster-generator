from ase import Atoms
from cluster_generator.filters.duplicates import get_metal_fingerprint, is_duplicate

def test_identical_structures_are_duplicates():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = atoms_1.copy()
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is True


def test_translated_structures_are_duplicates():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[4.0, 3.0, 1.0], [6.5, 3.0, 1.0], [5.25, 5.0, 1.0]])
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is True


def test_different_geometry_is_not_duplicate():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 3.0, 0.0]])
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is False


def test_fingerprint_ignores_rotation():
    atoms_1 = Atoms(symbols=["Co", "Co"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co"], positions=[[0.0, 0.0, 0.0], [0.0, 2.5, 0.0]])
    assert get_metal_fingerprint(atoms_1, {"Co"}) == get_metal_fingerprint(atoms_2, {"Co"})