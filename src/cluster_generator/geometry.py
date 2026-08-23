import numpy as np
from ase.data import atomic_numbers, covalent_radii


def get_preferred_distance(element_i, element_j, preferred_distances):
    '''
    Return the preferred starting distance for placing a new metal atom.
    '''
    pair_type = tuple(sorted((element_i, element_j)))

    if pair_type in preferred_distances:
        return preferred_distances[pair_type]

    radius_i = covalent_radii[atomic_numbers[element_i]]
    radius_j = covalent_radii[atomic_numbers[element_j]]
    return radius_i + radius_j


def get_growth_directions(atoms, n_directions=8, geometry="surface"):
    '''
    Generate surface-parallel, 3D or combined growth directions.
    '''
    if geometry == "both":
        return get_growth_directions(atoms, n_directions, "surface") + get_growth_directions(atoms, n_directions, "3d")

    if geometry == "surface":
        basis_1 = np.array(atoms.cell[0], dtype=float)
        basis_2 = np.array(atoms.cell[1], dtype=float)

        if np.linalg.norm(basis_1) == 0 or np.linalg.norm(basis_2) == 0:
            basis_1 = np.array([1.0, 0.0, 0.0])
            basis_2 = np.array([0.0, 1.0, 0.0])

        basis_1 = basis_1 / np.linalg.norm(basis_1)
        basis_2 = basis_2 - np.dot(basis_2, basis_1) * basis_1
        basis_2 = basis_2 / np.linalg.norm(basis_2)

        angles = np.linspace(0, 2 * np.pi, n_directions, endpoint=False)
        return [np.cos(angle) * basis_1 + np.sin(angle) * basis_2 for angle in angles]

    if geometry == "3d":
        golden_angle = np.pi * (3 - np.sqrt(5))
        directions = []

        for index in range(n_directions):
            z = 1 - 2 * (index + 0.5) / n_directions
            radius = np.sqrt(max(0.0, 1 - z * z))
            angle = golden_angle * index
            directions.append(np.array([radius * np.cos(angle), radius * np.sin(angle), z]))

        return directions

    raise ValueError("geometry must be 'surface', '3d' or 'both'")
