import argparse
from .io import get_input_files
from .workflows.growth import run_topological_growth
from pathlib import Path
from .io import get_input_files, read_structure, write_structure
from .workflows.growth import run_topological_growth, run_recursive_growth



def build_parser():
    '''
    Build the command-line interface.
    '''
    parser = argparse.ArgumentParser(prog="cluster-gen", description="Generate supported mono- and multimetallic cluster structures.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    grow = subparsers.add_parser("grow", help="Generate topological cluster-growth candidates.")
    grow.add_argument("inputs", nargs="+", help="Input files, directories or glob patterns.")
    grow.add_argument("--metals", nargs="+", required=True, help="Elements treated as the existing metal cluster.")
    grow.add_argument("--add", nargs="+", required=True, help="Element(s) allowed to be added during growth.")
    grow.add_argument("--geometry", choices=["surface", "3d", "both"], default="surface", help="Candidate direction mode.")
    grow.add_argument("--directions", type=int, default=8, help="Number of directions sampled per growth centre.")
    grow.add_argument("--target-size", type=int, help="Recursively grow clusters to this total number of metal atoms.")
    grow.add_argument("--size-range", nargs=2, type=int, metavar=("MIN", "MAX"), help="Generate and save all cluster sizes in this range.")
    grow.add_argument("--output-dir", default="generated_samples", help="Directory for generated structures.")

    return parser


def main():
    '''
    Run the command-line application.
    '''
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "grow":
        input_files = get_input_files(args.inputs)

        if not input_files:
            raise SystemExit("ERROR: No supported input files were found.")

        metals = set(args.metals)
        total_outputs = 0

        for input_file in input_files:
            try:
                if args.size_range:
                    min_size, max_size = args.size_range

                    atoms = read_structure(input_file)
                    structures, generations = run_recursive_growth(atoms, metals, args.add, max_size, args.geometry, args.directions, keep_generations=True)

                    output_dir = Path(args.output_dir)
                    output_dir.mkdir(parents=True, exist_ok=True)

                    for size in range(min_size, max_size + 1):
                        if size not in generations:
                            continue

                        for i, structure in enumerate(generations[size]):
                            output_name = output_dir / f"{input_file.stem}_n{size}_{i:04d}.extxyz"
                            write_structure(output_name, structure)
                            total_outputs += 1

                        print(f"{input_file.name}: saved {len(generations[size])} structure(s) at size {size}")
                elif args.target_size:
                    atoms = read_structure(input_file)
                    structures = run_recursive_growth(atoms, metals, args.add, args.target_size, args.geometry, args.directions)

                    output_dir = Path(args.output_dir)
                    output_dir.mkdir(parents=True, exist_ok=True)

                    for i, structure in enumerate(structures):
                        output_name = output_dir / f"{input_file.stem}_n{args.target_size}_{i:04d}.extxyz"
                        write_structure(output_name, structure)
                        total_outputs += 1

                    print(f"{input_file.name}: generated {len(structures)} structure(s) at size {args.target_size}")

                else:
                    for new_element in args.add:
                        structure_data, outputs = run_topological_growth(input_file, metals, new_element, args.geometry, args.directions, args.output_dir)
                        total_outputs += len(outputs)

                        print(f"{input_file.name}: add {new_element}, {structure_data['n_metals']} metal atom(s), "
                              f"{len(structure_data['growth_centres'])} centre(s), {len(structure_data['growth_edges'])} edge(s), "
                              f"{len(structure_data['growth_triangles'])} triangle(s), {len(outputs)} candidate structure(s)")

            except Exception as error:
                print(f"WARNING: Could not process {input_file}: {error}")

        print(f"Generated {total_outputs} candidate structure(s).")


if __name__ == "__main__":
    main()
