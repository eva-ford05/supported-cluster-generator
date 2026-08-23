from ase import Atom
from ..geometry import get_growth_directions
from ..distances import get_preferred_distance
from ..surface import get_surface_normal
import numpy as np

def generate_centre_growth(atoms, centre_index, new_element, preferred_distances, n_directions=8, geometry="surface"):
    '''
    Generate candidate positions for adding one metal atom around one existing metal centre.
    '''
    centre_element = atoms[centre_index].symbol
    centre_position = atoms[centre_index].position
    target_distance = get_preferred_distance(centre_element, new_element, preferred_distances)
    directions = get_growth_directions(atoms, n_directions=n_directions, geometry=geometry)

    candidates = []

    for direction in directions:
        position = centre_position + target_distance * direction
        candidates.append({
            "growth_type": "centre",
            "geometry": geometry,
            "centre": centre_index,
            "new_element": new_element,
            "position": position,
            "target_distance": target_distance,
        })

    return candidates

def generate_edge_growth(atoms, edge, new_element, preferred_distances, geometry="surface"):
    '''
    Generate candidate positions for adding one metal atom around an existing metal-metal edge
    '''
    i, j = edge

    position_1 = atoms[i].position
    position_2 = atoms[j].position

    element_1 = atoms[i].symbol
    element_2 = atoms[j].symbol

    distance_1 = get_preferred_distance(element_1, new_element, preferred_distances)
    distance_2 = get_preferred_distance(element_2, new_element, preferred_distances)

    edge_vector = position_2 - position_1
    edge_length = np.linalg.norm(edge_vector)

    if edge_length == 0:
        return []

    edge_unit = edge_vector / edge_length

    x = (distance_1**2 - distance_2**2 + edge_length**2) / (2 * edge_length)

    height_squared = distance_1**2 - x**2

    if height_squared < 0:
        return []

    height = np.sqrt(height_squared)
    circle_centre = position_1 + x * edge_unit

    candidates = []

    if geometry == "surface":
        surface_normal = get_surface_normal(atoms)

        perpendicular = np.cross(surface_normal, edge_unit)

        if np.linalg.norm(perpendicular) == 0:
            return []

        perpendicular = perpendicular / np.linalg.norm(perpendicular)

        positions = [
            circle_centre + height * perpendicular,
            circle_centre - height * perpendicular
        ]
    elif geometry == "3d":
        surface_normal = get_surface_normal(atoms)

        in_plane = np.cross(surface_normal, edge_unit)

        if np.linalg.norm(in_plane) == 0:
            return []
        in_plane = in_plane / np.linalg.norm(in_plane)

        positions = [
            circle_centre + height * in_plane,
            circle_centre - height * in_plane,
            circle_centre + height * surface_normal,
            circle_centre - height * surface_normal
        ]

    else:
        raise ValueError("geometry must be 'surface' or '3d'")

    for position in positions:
        candidates.append({"growth_type": "edge", "edge": edge, "new_element": new_element,
        "position": position, "distance_1": distance_1, "distance_2": distance_2})

    return candidates   

def generate_face_growth(atoms, triangle, new_element, preferred_distances):
    '''
    Generate candidate positions above and below an existing metal triangle
    '''

    i, j, k = triangle

    p1 = atoms[i].position
    p2 = atoms[j].position
    p3 = atoms[k].position

    v1 = p2 - p1
    v2 = p3 - p1

    normal = np.cross(v1, v2)

    if np.linalg.norm(normal) == 0:
        return []

    normal = normal / np.linalg.norm(normal)

    triangle_centre = (p1 + p2 + p3) / 3

    d1 = get_preferred_distance(atoms[i].symbol, new_element, preferred_distances)
    d2 = get_preferred_distance(atoms[j].symbol, new_element, preferred_distances)
    d3 = get_preferred_distance(atoms[k].symbol, new_element, preferred_distances)

    target_distance = (d1 + d2 + d3) / 3

    side_1 = np.linalg.norm(p2 - p1)
    side_2 = np.linalg.norm(p3 - p2)
    side_3 = np.linalg.norm(p1 - p3)

    average_side = (side_1 + side_2 + side_3) / 3
    in_plane_radius = average_side / np.sqrt(3)

    height_squared = target_distance**2 - in_plane_radius**2

    if height_squared < 0:
        return []

    height = np.sqrt(height_squared)

    positions = [
        triangle_centre + height * normal,
        triangle_centre - height * normal
    ]

    candidates = []

    for position in positions: 
        candidates.append({"growth_type": "face", "triangle": triangle, "new_element": new_element,
        "position": position, "target_distance": target_distance})

    return candidates

def build_candidate_structure(atoms, candidate):
    '''
    Copy the parent structure and append the candidate metal atom.
    '''
    new_atoms = atoms.copy()
    new_atoms.append(Atom(candidate["new_element"], position=candidate["position"]))
    return new_atoms
