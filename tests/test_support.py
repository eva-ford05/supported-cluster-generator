from ase import Atoms
from cluster_generator.filters.support import is_inside_support


def test_atom_above_support_is_allowed():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 7.0]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}) is False


def test_atom_slightly_below_surface_is_allowed():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.7]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is False


def test_atom_inside_support_is_rejected():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.0]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is True


def test_tolerance_changes_penetration_cutoff():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.7]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is False
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.2) is True