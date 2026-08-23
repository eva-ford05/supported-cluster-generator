import numpy as np
from ase import Atoms
from cluster_generator.generators.topological import generate_edge_growth

def test_edge_growth_surface_positions():
    atoms = Atoms(
        symbols=["Co", "Co"],
        positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=True,)

    preferred_distances = {("Co", "Ru"): 2.6}

    candidates = generate_edge_growth(
        atoms, (0, 1), "Ru", preferred_distances, geometry="surface")

    assert len(candidates) == 2

    for candidate in candidates:
        position = candidate["position"]

        distance_1 = np.linalg.norm(position - atoms[0].position)
        distance_2 = np.linalg.norm(position - atoms[1].position)

        assert np.isclose(distance_1, 2.6)
        assert np.isclose(distance_2, 2.6)
        assert np.isclose(position[2], 0.0)