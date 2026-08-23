default_metals = {"Co", "Mn", "Ru"}
default_extensions = {".xyz", ".extxyz", ".in"}

# Temporary graph-neighbour cutoffs. These can be replaced with better chemistry-specific values later.
bond_cutoffs = {
    ("Co", "Co"): 3.6,
    ("Co", "Mn"): 3.6,
    ("Co", "Ru"): 3.6,
    ("Mn", "Mn"): 3.6,
    ("Mn", "Ru"): 3.6,
    ("Ru", "Ru"): 3.6,
}

bond_tolerance = 1.20

# Optional placement-distance overrides. ASE covalent radii are used when a pair is not listed.
preferred_distances = {}
