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


def is_duplicate_fingerprint(atoms, seen_fingerprints, metals, decimals=3):
    '''
    Return True if the metal fingerprint has already been seen.
    '''

    fingerprint = get_metal_fingerprint(atoms, metals, decimals)

    if fingerprint in seen_fingerprints:
        return True

    seen_fingerprints.add(fingerprint)
    return False

def get_distance_matrix_fingerprint(atoms, metals, decimals=3):
    '''
    Return an element aware metal distance matrix fingerprint
    '''

    metal_indices = [atom.index for atom in atoms if atom.symbol in metals]
    rows = []

    for i in metal_indices:
        distances = []

        for j in metal_indices:
            if i == j:
                continue
            pair = tuple(sorted((atoms[i].symbol, atoms[j].symbol)))

            distance = round(float(atoms.get_distance(i, j, mic=True)), decimals)
            distances.append((pair, distance))

        rows.append((atoms[i].symbol, tuple(sorted(distances))))

    return tuple(sorted(rows))

def get_support_environment_fingerprint(atoms, metals, cutoff=3.0, decimals=3):
    '''
    Return a fingerprint describing the local support environment around the metal cluster.
    '''

    metal_indices = [atom.index for atom in atoms if atom.symbol in metals]
    support_indices = [atom.index for atom in atoms if atom.symbol not in metals]

    environments = []

    for i in metal_indices:
        neighbours = []

        for j in support_indices:
            distance = float(atoms.get_distance(i, j, mic=True))

            if distance <= cutoff:
                neighbours.append((atoms[j].symbol, round(distance, decimals)))

        environments.append((atoms[i].symbol, tuple(sorted(neighbours))))

    return tuple(sorted(environments))

def get_supported_cluster_fingerprint(atoms, metals, support_cutoff=3.0, decimals=3):
    '''
    Return a fingerprint containing both cluster geometry and local support environment.
    '''

    metal_fingerprint = get_distance_matrix_fingerprint(atoms, metals, decimals=decimals)
    support_fingerprint = get_support_environment_fingerprint(atoms, metals, cutoff=support_cutoff, decimals=decimals)

    return metal_fingerprint, support_fingerprint

def is_duplicate_supported(atoms, seen_fingerprints, metals, support_cutoff=3.0, decimals=3):
    '''
    Return True if the supported cluster fingerprint has already been seen.
    '''

    fingerprint = get_supported_cluster_fingerprint(atoms, metals, support_cutoff=support_cutoff, decimals=decimals)

    if fingerprint in seen_fingerprints:
        return True

    seen_fingerprints.add(fingerprint)
    return False