import numpy as np
from ase import Atoms
from cluster_generator.generators.topological import generate_face_growth

def test_face_growth_positions():
    side = 2.5
    height_2d = np.sqrt(3) * side / 2

    atoms = Atoms(symbols=["Co", "Co", "Co"], positions=[
            [0.0, 0.0, 0.0], [side, 0.0, 0.0], [side / 2, height_2d, 0.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=True)

    preferred_distances = {("Co", "Ru"): 2.6}

    candidates = generate_face_growth(atoms, (0, 1, 2), "Ru", preferred_distances)

    assert len(candidates) == 2

    assert candidates[0]["position"][2] > 0
    assert candidates[1]["position"][2] < 0