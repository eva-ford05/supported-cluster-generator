from ase import Atoms
from cluster_generator.filters.support import is_inside_support, get_local_surface_height, get_surface_normal


def test_atom_above_support_is_allowed():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 7.0]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}) is False


def test_atom_slightly_below_surface_is_allowed():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.7]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is False


def test_atom_inside_support_is_rejected():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.0]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is True


def test_tolerance_changes_penetration_cutoff():
    atoms = Atoms(symbols=["Ti", "O", "Ru"], positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [0.5, 0.5, 5.7]], cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]], pbc=True)
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5) is False
    assert is_inside_support(atoms, {"Ru"}, tolerance=0.2) is True

def test_local_surface_height_uses_nearby_atoms():
    atoms = Atoms(
        symbols=["O", "Ti", "O"],
        positions=[[0.0, 0.0, 5.0], [1.0, 0.0, 6.0], [8.0, 8.0, 9.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=True,
    )

    height = get_local_surface_height(atoms, {"Ru"}, [0.5, 0.5, 7.0], radius=3.0)

    assert height == 6.0

def test_local_surface_height_changes_with_position():
    atoms = Atoms(
        symbols=["O", "O"],
        positions=[[1.0, 1.0, 5.0], [7.0, 7.0, 8.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=True,
    )

    low_region = get_local_surface_height(atoms, {"Ru"}, [1.0, 1.0, 7.0], radius=2.0)
    high_region = get_local_surface_height(atoms, {"Ru"}, [7.0, 7.0, 10.0], radius=2.0)

    assert low_region == 5.0
    assert high_region == 8.0

def test_distant_high_support_atom_does_not_reject_candidate():
    atoms = Atoms(
        symbols=["O", "Ti", "O", "Ru"],
        positions=[
            [0.0, 0.0, 5.0],
            [1.0, 0.0, 6.0],
            [8.0, 8.0, 9.0],
            [0.5, 0.5, 6.2],
        ],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=True,
    )

    assert is_inside_support(atoms, {"Ru"}, tolerance=0.5, local_radius=3.0) is False

def test_local_surface_height_respects_periodic_boundaries():
    atoms = Atoms(
        symbols=["O", "O"],
        positions=[[0.2, 5.0, 6.0], [5.0, 5.0, 8.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=[True, True, False],
    )

    height = get_local_surface_height(atoms, {"Ru"}, [9.8, 5.0, 7.0], radius=1.0)

    assert height == 6.0

def test_surface_normal_points_towards_cluster():
    atoms = Atoms(
        symbols=["Ti", "O", "Ru"],
        positions=[[0.0, 0.0, 4.0], [1.0, 0.0, 5.0], [0.5, 0.5, 7.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=[True, True, False],
    )

    normal = get_surface_normal(atoms, {"Ru"})

    assert normal[2] > 0

def test_surface_normal_flips_for_cluster_below_support():
    atoms = Atoms(
        symbols=["Ti", "O", "Ru"],
        positions=[[0.0, 0.0, 6.0], [1.0, 0.0, 5.0], [0.5, 0.5, 3.0]],
        cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]],
        pbc=[True, True, False],
    )

    normal = get_surface_normal(atoms, {"Ru"})

    assert normal[2] < 0