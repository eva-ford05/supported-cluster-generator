import numpy as np

def get_surface_normal(atoms):
    '''
    Return a unit vector normal to the support surface.
    '''

    normal = np.cross(atoms.cell[0], atoms.cell[1])
    return normal / np.linalg.norm(normal)

def get_surface_basis(atoms):
    '''
    Return two orthonormal vectors spanning the support surface.
    '''

    basis_1 = np.array(atoms.cell[0], dtype=float)
    basis_2 = np.array(atoms.cell[1], dtype=float)

    basis_1 = basis_1 / np.linalg.norm(basis_1)

    basis_2 = basis_2 - np.dot(basis_2, basis_1) * basis_1
    basis_2 = basis_2 / np.linalg.norm(basis_2)

    return basis_1, basis_2