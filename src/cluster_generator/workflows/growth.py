from pathlib import Path
from ..config import bond_cutoffs, bond_tolerance, preferred_distances
from ..generators.topological import build_candidate_structure, generate_centre_growth, generate_edge_growth, generate_face_growth
from ..io import read_structure, write_structure
from ..topology import analyse_structure


def run_topological_growth(input_file, metals, new_element, geometry="surface", n_directions=8, output_dir="generated_samples"):
    '''
    Run one-centre topological growth for one input structure.
    '''
    input_file = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = read_structure(input_file)
    structure_data = analyse_structure(atoms, metals, bond_cutoffs, bond_tolerance)
    outputs = []

    # 1-centre growth
    for centre in structure_data["growth_centres"]:
        candidates = generate_centre_growth(atoms, centre, new_element, preferred_distances, n_directions=n_directions, geometry=geometry)

        for candidate_number, candidate in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidate)
            output_name = output_dir / f"{input_file.stem}_centre{centre}_{new_element}_{geometry}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)

    # 2-centre growth
    for edge in structure_data["growth_edges"]:
        candidates = generate_edge_growth(atoms, edge, new_element, preferred_distances, geometry=geometry)
        print(f"Edge {edge}: {len(candidates)} candidate(s)")
        
        for candidate_number, candidates in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidates)

            output_name = output_dir / f"{input_file.stem}_edge{edge[0]}-{edge[1]}_{new_element}_{geometry}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)

    # 3-centre growth
    for triangle in structure_data["growth_triangles"]:
        candidates = generate_face_growth(atoms, triangle, new_element, preferred_distances)

        for candidate_number, candidate in enumerate(candidates):
            new_atoms = build_candidate_structure(atoms, candidate)
            output_name = output_dir / f"{input_file.stem}_face{triangle[0]}-{triangle[1]}-{triangle[2]}_{new_element}_{candidate_number}.extxyz"
            write_structure(output_name, new_atoms)
            outputs.append(output_name)
            
    return structure_data, outputs
