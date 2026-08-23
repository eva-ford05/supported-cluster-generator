from ase import Atom
from ..geometry import get_growth_directions, get_preferred_distance


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


def build_candidate_structure(atoms, candidate):
    '''
    Copy the parent structure and append the candidate metal atom.
    '''
    new_atoms = atoms.copy()
    new_atoms.append(Atom(candidate["new_element"], position=candidate["position"]))
    return new_atoms
