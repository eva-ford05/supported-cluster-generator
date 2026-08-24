from pathlib import Path
from ..config import bond_cutoffs, bond_tolerance, preferred_distances, metal_clash_scale, support_clash_scale, support_local_radius, support_penetration_tolerance
from ..generators.topological import build_candidate_structure, generate_centre_growth, generate_edge_growth, generate_face_growth
from ..io import read_structure, write_structure
from ..topology import analyse_structure
from ..filters.clashes import has_atomic_clash
from ..filters.support import is_inside_support
from ..filters.duplicates import is_duplicate_fingerprint

def generate_topological_children(atoms, metals, new_elements, geometry="surface", n_directions=8):
    '''
    Generate all valid N -> N+1 child structures for the allowed added elements.
    '''

    structure_data = analyse_structure(atoms, metals, bond_cutoffs, bond_tolerance)

    children = []
    seen_fingerprints = set()
    stats = {"raw": 0, "clash": 0, "support": 0, "within_parent_duplicate": 0, "kept": 0}

    for new_element in new_elements:

        # 1-centre growth
        for centre in structure_data["growth_centres"]:
            candidates = generate_centre_growth(atoms, centre, new_element, preferred_distances, n_directions=n_directions, geometry=geometry)

            for candidate in candidates:
                stats["raw"] += 1
                new_atoms = build_candidate_structure(atoms, candidate)

                if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                    stats["clash"] += 1
                    continue

                if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                    stats["support"] += 1
                    continue

                if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                    stats["within_parent_duplicate"] += 1
                    continue

                children.append(new_atoms)
                stats["kept"] += 1

        # 2-centre growth
        for edge in structure_data["growth_edges"]:
            candidates = generate_edge_growth(atoms, edge, new_element, preferred_distances, geometry=geometry)

            for candidate in candidates:
                stats["raw"] += 1
                new_atoms = build_candidate_structure(atoms, candidate)

                if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                    stats["clash"] += 1
                    continue

                if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                    stats["support"] += 1
                    continue

                if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                    stats["within_parent_duplicate"] += 1
                    continue

                children.append(new_atoms)
                stats["kept"] += 1

        # 3-centre growth
        for triangle in structure_data["growth_triangles"]:
            candidates = generate_face_growth(atoms, triangle, new_element, preferred_distances)

            for candidate in candidates:
                stats["raw"] += 1
                new_atoms = build_candidate_structure(atoms, candidate)

                if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                    stats["clash"] += 1
                    continue

                if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                    stats["support"] += 1
                    continue

                if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                    stats["within_parent_duplicate"] += 1
                    continue

                children.append(new_atoms)
                stats["kept"] += 1

    return children, stats


def run_recursive_growth(atoms, metals, new_elements, target_size, geometry="surface", n_directions=8, keep_generations=False):
    '''
    Recursively grow a supported cluster until the target metal count is reached.
    '''

    current_structures = [atoms]
    current_size = len([atom for atom in atoms if atom.symbol in metals])
    generations = {}

    while current_size < target_size:
        next_structures = []
        seen_fingerprints = set()

        generation_stats = {
            "raw": 0,
            "clash": 0,
            "support": 0,
            "within_parent_duplicate": 0,
            "cross_parent_duplicate": 0,
        }

        for parent in current_structures:
            children, stats = generate_topological_children(parent, metals, new_elements, geometry=geometry, n_directions=n_directions)

            generation_stats["raw"] += stats["raw"]
            generation_stats["clash"] += stats["clash"]
            generation_stats["support"] += stats["support"]
            generation_stats["within_parent_duplicate"] += stats["within_parent_duplicate"]

            for child in children:
                if is_duplicate_fingerprint(child, seen_fingerprints, metals):
                    generation_stats["cross_parent_duplicate"] += 1
                    continue

                next_structures.append(child)

        if not next_structures:
            break

        current_structures = next_structures
        current_size += 1

        if keep_generations:
            generations[current_size] = current_structures.copy()

        print(f"size {current_size}: {generation_stats['raw']} raw, {generation_stats['clash']} clash, "
              f"{generation_stats['support']} support, {generation_stats['within_parent_duplicate']} within-parent duplicate, "
              f"{generation_stats['cross_parent_duplicate']} cross-parent duplicate, {len(current_structures)} kept")

    if keep_generations:
        return current_structures, generations

    return current_structures


def run_topological_growth(input_file, metals, new_element, geometry="surface", n_directions=8, output_dir="generated_samples"):
    '''
    Run topological cluster growth for one input structure.
    '''

    input_file = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = read_structure(input_file)
    structure_data = analyse_structure(atoms, metals, bond_cutoffs, bond_tolerance)

    outputs = []
    seen_fingerprints = set()

    # 1-centre growth
    for centre in structure_data["growth_centres"]:
        candidates = generate_centre_growth(atoms, centre, new_element, preferred_distances, n_directions=n_directions, geometry=geometry)

        for candidate_number, candidate in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidate)

            if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                continue

            if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                continue

            if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                continue

            output_name = output_dir / f"{input_file.stem}_centre{centre}_{new_element}_{geometry}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)

    # 2-centre growth
    for edge in structure_data["growth_edges"]:
        candidates = generate_edge_growth(atoms, edge, new_element, preferred_distances, geometry=geometry)

        for candidate_number, candidate in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidate)

            if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                continue

            if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                continue

            if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                continue

            output_name = output_dir / f"{input_file.stem}_edge{edge[0]}-{edge[1]}_{new_element}_{geometry}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)

    # 3-centre growth
    for triangle in structure_data["growth_triangles"]:
        candidates = generate_face_growth(atoms, triangle, new_element, preferred_distances)

        for candidate_number, candidate in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidate)

            if has_atomic_clash(new_atoms, metals, metal_scale=metal_clash_scale, support_scale=support_clash_scale):
                continue

            if is_inside_support(new_atoms, metals, tolerance=support_penetration_tolerance, local_radius=support_local_radius):
                continue

            if is_duplicate_fingerprint(new_atoms, seen_fingerprints, metals):
                continue

            output_name = output_dir / f"{input_file.stem}_face{triangle[0]}-{triangle[1]}-{triangle[2]}_{new_element}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)

    return structure_data, outputs