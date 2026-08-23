from ..distances import get_clash_distance

def has_atomic_clash(atoms, metals, new_atom_index=-1, clash_overrides=None, metal_scale=0.60, support_scale=0.50):
    '''
    Return True only for severe atomic overlaps.

    Minimum distances are estimated from pair-specific covalent radii,
    with separate scaling for metal-metal and metal-support contacts.
    '''

    if new_atom_index < 0:
        new_atom_index = len(atoms) + new_atom_index

    new_atom = atoms[new_atom_index]

    for atom in atoms:
        if atom.index == new_atom_index:
            continue

        distance = atoms.get_distance(new_atom_index, atom.index, mic=True)

        minimum_distance = get_clash_distance(new_atom.symbol, atom.symbol, metals, clash_overrides, metal_scale, support_scale)

        if distance < minimum_distance:
            return True

    return False