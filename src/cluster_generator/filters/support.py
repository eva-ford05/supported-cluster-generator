import numpy as np

def is_inside_support(atoms, metals, new_atom_index=-1, tolerance=0.5):
    '''
    Return True if the newly added metal sits too far inside the support

    The support surface is estimated from the highest support atoms along the
    slab normal
    '''

    if new_atom_index < 0:
        new_atom_index = len(atoms) + new_atom_index

    surface_normal = np.cross(atoms.cell[0], atoms.cell[1])
    surface_normal = surface_normal / np.linalg.norm(surface_normal)

    support_indices = [atom.index for atom in atoms if atom.symbol not in metals and atom.index 
    != new_atom_index]

    if not support_indices:
        return False

    support_heights = [np.dot(atoms[i].position, surface_normal) for i in support_indices]
    surface_height = max(support_heights)

    new_height = np.dot(atoms[new_atom_index].position, surface_normal)

    return bool(new_height < surface_height - tolerance)