from itertools import combinations

def get_metal_fingerprint(atoms, metals, decimals=3):
    '''
    Return a simple metal-cluster fingerprint based on element-labelled pair distances.
    '''

    metal_indices = [atom.index for atom in atoms if atom.symbol in metals]
    fingerprint = []

    for i, j in combinations(metal_indices, 2):
        pair = tuple(sorted((atoms[i].symbol, atoms[j].symbol)))
        distance = round(float(atoms.get_distance(i, j, mic=True)), decimals)
        fingerprint.append((pair, distance))

    return tuple(sorted(fingerprint))


def is_duplicate(atoms, existing_structures, metals, decimals=3):
    '''
    Return True if the metal framework matches one already kept.
    '''

    fingerprint = get_metal_fingerprint(atoms, metals, decimals)

    for existing_atoms in existing_structures:
        existing_fingerprint = get_metal_fingerprint(existing_atoms, metals, decimals)

        if fingerprint == existing_fingerprint:
            return True

    return False