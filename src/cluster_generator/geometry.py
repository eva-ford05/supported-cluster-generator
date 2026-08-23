import numpy as np
from .surface import get_surface_basis

def get_growth_directions(atoms, n_directions=8, geometry="surface"):
    '''
    Generate surface-parallel, 3D or combined growth directions.
    '''

    if geometry == "both":
        return get_growth_directions(atoms, n_directions, "surface") + get_growth_directions(atoms, n_directions, "3d")

    if geometry == "surface":
        basis_1, basis_2 = get_surface_basis(atoms)

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
