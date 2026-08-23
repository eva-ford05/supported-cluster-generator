from ase.data import atomic_numbers, covalent_radii

def get_pair_radius_sum(element_i, element_j):
    '''
    Return the usm of the ASE covalent radii for two elements
    '''

    radius_i = covalent_radii[atomic_numbers[element_i]]
    radius_j = covalent_radii[atomic_numbers[element_j]]

    return radius_i + radius_j

def get_bond_cutoff(element_i, element_j, bond_cutoffs=None, bond_tolerance=1.20):
    '''
    Return the cutoff used to decide whether two metals are connected.

    A user-defined pair cutoff is used when available. Otherwise the cutoff
    is estimated from the sum of ASE covalent radii.
    '''

    pair_type = tuple(sorted((element_i, element_j)))

    if bond_cutoffs and pair_type in bond_cutoffs:
        return bond_cutoffs[pair_type]

    return bond_tolerance * get_pair_radius_sum(element_i, element_j)

def get_preferred_distance(element_i, element_j, preferred_distances=None):
    '''
    Return the preferred starting distance used when placing a new metal atom.

    A user-defined pair distance is used when available. Otherwise the sum of ASE covalent radii is used.
    '''

    pair_type = tuple(sorted((element_i, element_j)))

    if preferred_distances and pair_type in preferred_distances:
        return preferred_distances[pair_type]

    return get_pair_radius_sum(element_i, element_j)

def get_clash_distance(element_i, element_j, metals, clash_overrides=None, metal_scale=0.60, support_scale=0.50):
    '''
    Return the minimum allowed distance before two atoms are treated as a severe clash.

    Metal-metal and metal-support pairs use different scaling factors.
    Pair-specific overrides can be supplied when needed.
    '''

    pair_type = tuple(sorted((element_i, element_j)))

    if clash_overrides and pair_type in clash_overrides:
        return clash_overrides[pair_type]

    radius_sum = get_pair_radius_sum(element_i, element_j)

    if element_j in metals:
        return metal_scale * radius_sum

    return support_scale * radius_sum