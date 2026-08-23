import numpy as np
from cluster_generator.geometry import get_growth_directions


class DummyAtoms:
    def __init__(self):
        self.cell = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]])


def test_surface_directions_have_zero_z():
    atoms = DummyAtoms()
    directions = get_growth_directions(atoms, n_directions=8, geometry="surface")

    assert len(directions) == 8

    for direction in directions:
        assert np.isclose(direction[2], 0.0)


def test_3d_directions_are_unit_vectors():
    atoms = DummyAtoms()
    directions = get_growth_directions(atoms, n_directions=16, geometry="3d")

    assert len(directions) == 16

    for direction in directions:
        assert np.isclose(np.linalg.norm(direction), 1.0)
