from ase import Atoms
from cluster_generator.filters.duplicates import get_metal_fingerprint, is_duplicate, is_duplicate_fingerprint, get_distance_matrix_fingerprint, get_support_environment_fingerprint, get_supported_cluster_fingerprint

def test_identical_structures_are_duplicates():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = atoms_1.copy()
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is True


def test_translated_structures_are_duplicates():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[4.0, 3.0, 1.0], [6.5, 3.0, 1.0], [5.25, 5.0, 1.0]])
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is True


def test_different_geometry_is_not_duplicate():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 3.0, 0.0]])
    assert is_duplicate(atoms_2, [atoms_1], {"Co", "Ru"}) is False


def test_fingerprint_ignores_rotation():
    atoms_1 = Atoms(symbols=["Co", "Co"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co"], positions=[[0.0, 0.0, 0.0], [0.0, 2.5, 0.0]])
    assert get_metal_fingerprint(atoms_1, {"Co"}) == get_metal_fingerprint(atoms_2, {"Co"})

def test_duplicate_fingerprint_is_added_and_detected():
    atoms = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])

    seen_fingerprints = set()

    assert is_duplicate_fingerprint(atoms, seen_fingerprints, {"Co", "Ru"}) is False
    assert is_duplicate_fingerprint(atoms, seen_fingerprints, {"Co", "Ru"}) is True

def test_different_fingerprints_are_kept():
    atoms_1 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Co", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 3.0, 0.0]])

    seen_fingerprints = set()

    assert is_duplicate_fingerprint(atoms_1, seen_fingerprints, {"Co", "Ru"}) is False
    assert is_duplicate_fingerprint(atoms_2, seen_fingerprints, {"Co", "Ru"}) is False

def test_distance_matrix_fingerprint_is_translation_invariant():
    atoms_1 = Atoms(symbols=["Co", "Mn", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.0, 2.0, 0.0]])
    atoms_2 = atoms_1.copy()
    atoms_2.translate([5.0, 4.0, 3.0])

    assert get_distance_matrix_fingerprint(atoms_1, {"Co", "Mn", "Ru"}) == get_distance_matrix_fingerprint(atoms_2, {"Co", "Mn", "Ru"})

def test_distance_matrix_fingerprint_detects_different_geometry():
    atoms_1 = Atoms(symbols=["Co", "Mn", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.0, 2.0, 0.0]])
    atoms_2 = Atoms(symbols=["Co", "Mn", "Ru"], positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.0, 3.0, 0.0]])

    assert get_distance_matrix_fingerprint(atoms_1, {"Co", "Mn", "Ru"}) != get_distance_matrix_fingerprint(atoms_2, {"Co", "Mn", "Ru"})

def test_support_environment_fingerprint_detects_different_sites():
    atoms_1 = Atoms(symbols=["O", "Ti", "Ru"], positions=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 2.0]])

    atoms_2 = Atoms(symbols=["O", "Ti", "Ru"], positions=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.5, 0.0, 2.0]])

    assert get_support_environment_fingerprint(atoms_1, {"Ru"}) != get_support_environment_fingerprint(atoms_2, {"Ru"})

def test_supported_fingerprint_matches_identical_structure():
    atoms_1 = Atoms(symbols=["O", "Ti", "Co", "Ru"],
        positions=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0]])

    atoms_2 = atoms_1.copy()

    assert get_supported_cluster_fingerprint(atoms_1, {"Co", "Ru"}) == get_supported_cluster_fingerprint(atoms_2, {"Co", "Ru"})

def test_supported_fingerprint_distinguishes_support_environment():
    atoms_1 = Atoms(symbols=["O", "Ti", "Co", "Ru"],
        positions=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0]])

    atoms_2 = Atoms(symbols=["O", "Ti", "Co", "Ru"],
        positions=[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.5, 0.0, 2.0], [2.5, 0.0, 2.0]])

    assert get_supported_cluster_fingerprint(atoms_1, {"Co", "Ru"}) != get_supported_cluster_fingerprint(atoms_2, {"Co", "Ru"})