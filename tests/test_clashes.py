from ase import Atoms
from cluster_generator.filters.clashes import has_atomic_clash


def test_metal_metal_clash():
    atoms = Atoms(symbols=["Co", "Ru"], positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert has_atomic_clash(atoms, {"Co", "Ru", "Mn"}) is True


def test_reasonable_metal_metal_distance():
    atoms = Atoms(symbols=["Co", "Ru"], positions=[[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]])
    assert has_atomic_clash(atoms, {"Co", "Ru", "Mn"}) is False


def test_metal_support_clash():
    atoms = Atoms(symbols=["O", "Ru"], positions=[[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]])
    assert has_atomic_clash(atoms, {"Co", "Ru", "Mn"}) is True


def test_reasonable_metal_support_distance():
    atoms = Atoms(symbols=["O", "Ru"], positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    assert has_atomic_clash(atoms, {"Co", "Ru", "Mn"}) is False