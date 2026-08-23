from cluster_generator.distances import get_pair_radius_sum, get_bond_cutoff, get_preferred_distance, get_clash_distance

def test_pair_radius_sum_is_positive():
    distance = get_pair_radius_sum("Co", "Ru")
    assert distance > 0


def test_pair_radius_sum_is_symmetric():
    assert get_pair_radius_sum("Co", "Ru") == get_pair_radius_sum("Ru", "Co")


def test_different_pairs_have_different_radius_sums():
    co_ru = get_pair_radius_sum("Co", "Ru")
    ru_o = get_pair_radius_sum("Ru", "O")

    assert co_ru != ru_o

def test_bond_cutoff_uses_override():
    bond_cutoffs = {("Co", "Ru"): 3.2}
    assert get_bond_cutoff("Co", "Ru", bond_cutoffs) == 3.2


def test_bond_cutoff_override_is_symmetric():
    bond_cutoffs = {("Co", "Ru"): 3.2}
    assert get_bond_cutoff("Ru", "Co", bond_cutoffs) == 3.2


def test_bond_cutoff_falls_back_to_radius_sum():
    expected = 1.20 * get_pair_radius_sum("Co", "Ru")
    result = get_bond_cutoff("Co", "Ru", bond_cutoffs={})

    assert result == expected

def test_preferred_distance_uses_override():
    preferred_distances = {("Co", "Ru"): 2.55}
    assert get_preferred_distance("Co", "Ru", preferred_distances) == 2.55

def test_preferred_distance_override_is_symmetric():
    preferred_distances = {("Co", "Ru"): 2.55}
    assert get_preferred_distance("Ru", "Co", preferred_distances) == 2.55

def test_preferred_distance_falls_back_to_radius_sum():
    expected = get_pair_radius_sum("Co", "Ru")
    result = get_preferred_distance("Co", "Ru", preferred_distances={})
    assert result == expected

def test_metal_clash_distance_uses_radius_sum():
    expected = 0.60 * get_pair_radius_sum("Co", "Ru")
    result = get_clash_distance("Co", "Ru", {"Co", "Ru"})
    assert result == expected

def test_support_clash_distance_uses_radius_sum():
    expected = 0.50 * get_pair_radius_sum("Ru", "O")
    result = get_clash_distance("Ru", "O", {"Co", "Ru"})
    assert result == expected

def test_different_support_atoms_have_different_clash_distances():
    ru_o = get_clash_distance("Ru", "O", {"Ru"})
    ru_ti = get_clash_distance("Ru", "Ti", {"Ru"})
    assert ru_o != ru_ti

def test_clash_distance_uses_override():
    clash_overrides = {("O", "Ru"): 1.1}
    assert get_clash_distance("Ru", "O", {"Ru"}, clash_overrides) == 1.1

def test_clash_override_is_symmetric():
    clash_overrides = {("O", "Ru"): 1.1}
    assert get_clash_distance("O", "Ru", {"Ru"}, clash_overrides) == 1.1