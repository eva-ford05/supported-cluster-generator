import numpy as np
from ase.geometry import find_mic
from ..surface import get_surface_normal as get_base_surface_normal
from ..selection import get_cluster_indices, get_support_indices

def get_local_surface_height(atoms, metals, position, radius=3.0, surface_normal=None):
    '''
    Estimate the local support height around a position.

    Only support atoms within the given lateral radius are considered.
    '''

    if surface_normal is None:
        surface_normal = get_surface_normal(atoms, metals)

    support_indices = get_support_indices(atoms, metals)
    local_heights = []

    for i in support_indices:
        displacement = atoms[i].position - position
        displacement, _ = find_mic(displacement, atoms.cell, atoms.pbc)

        normal_component = np.dot(displacement, surface_normal) * surface_normal
        lateral_component = displacement - normal_component
        lateral_distance = np.linalg.norm(lateral_component)

        if lateral_distance <= radius:
            height = np.dot(atoms[i].position, surface_normal)
            local_heights.append(height)

    if not local_heights:
        return None

    return max(local_heights)

def is_inside_support(atoms, metals, new_atom_index=-1, tolerance=0.5, local_radius=3.0):
    '''
    Return True if the newly added metal sits too far inside the local support surface.
    '''

    if new_atom_index < 0:
        new_atom_index = len(atoms) + new_atom_index

    surface_normal = get_surface_normal(atoms, metals, exclude_index=new_atom_index)

    new_position = atoms[new_atom_index].position
    new_height = np.dot(new_position, surface_normal)

    surface_height = get_local_surface_height(atoms, metals, new_position, radius=local_radius, surface_normal=surface_normal)

    if surface_height is None:
        return False

    return bool(new_height < surface_height - tolerance)

def get_surface_normal(atoms, metals, exclude_index=None):
    '''
    Return the slab normal oriented towards the existing supported metal cluster.
    '''

    surface_normal = get_base_surface_normal(atoms)

    metal_indices = [i for i in get_cluster_indices(atoms, metals) if i != exclude_index]

    support_indices = get_support_indices(atoms, metals)

    if not metal_indices or not support_indices:
        return surface_normal

    metal_heights = [np.dot(atoms[i].position, surface_normal) for i in metal_indices]
    support_heights = [np.dot(atoms[i].position, surface_normal) for i in support_indices]

    if np.mean(metal_heights) < np.mean(support_heights):
        surface_normal = -surface_normal

    return surface_normal