
def has_atomic_clash(atoms, metals, new_atom_index=-1, metal_min_distance=1.6, support_min_distance=1.2):
    '''
    Return True only for severe atomic overlaps.

    Metal-metal and metal-support contacts use separate minimum distances.
    '''

    if new_atom_index < 0:
        new_atom_index = len(atoms) + new_atom_index

    new_atom = atoms[new_atom_index]

    for atom in atoms:
        if atom.index == new_atom_index:
            continue

        distance = atoms.get_distance(new_atom_index, atom.index, mic=True)

        if atom.symbol in metals:
            minimum_distance = metal_min_distance
        else:
            minimum_distance = support_min_distance

        if distance < minimum_distance:
            return True

    return False