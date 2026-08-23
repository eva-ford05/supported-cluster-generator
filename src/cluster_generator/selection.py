def get_cluster_indices(atoms, metals):
    '''
    Return atom indices treated as part of the supported cluster.
    '''

    return [atom.index for atom in atoms if atom.symbol in metals]


def get_support_indices(atoms, metals):
    '''
    Return atom indices treated as part of the support.
    '''

    return [atom.index for atom in atoms if atom.symbol not in metals]