import numpy as np
from ase import Atoms
from cluster_generator.surface import get_surface_normal, get_surface_basis


def test_surface_normal_is_unit_length():
    atoms = Atoms(cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    normal = get_surface_normal(atoms)

    assert np.isclose(np.linalg.norm(normal), 1.0)


def test_surface_normal_is_perpendicular_to_surface():
    atoms = Atoms(cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    normal = get_surface_normal(atoms)

    assert np.isclose(np.dot(normal, atoms.cell[0]), 0.0)
    assert np.isclose(np.dot(normal, atoms.cell[1]), 0.0)


def test_surface_basis_is_orthonormal():
    atoms = Atoms(cell=[[10.0, 0.0, 0.0], [2.0, 9.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    basis_1, basis_2 = get_surface_basis(atoms)

    assert np.isclose(np.linalg.norm(basis_1), 1.0)
    assert np.isclose(np.linalg.norm(basis_2), 1.0)
    assert np.isclose(np.dot(basis_1, basis_2), 0.0)