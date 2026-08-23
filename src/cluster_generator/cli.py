import argparse
from .config import default_metals
from .io import get_input_files
from .workflows.growth import run_topological_growth


def build_parser():
    '''
    Build the command-line interface.
    '''
    parser = argparse.ArgumentParser(prog="cluster-gen", description="Generate supported mono- and multimetallic cluster structures.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    grow = subparsers.add_parser("grow", help="Generate topological cluster-growth candidates.")
    grow.add_argument("inputs", nargs="+", help="Input files, directories or glob patterns.")
    grow.add_argument("--metals", nargs="+", default=sorted(default_metals), help="Elements treated as the existing metal cluster.")
    grow.add_argument("--add", required=True, help="Element to add during this growth step.")
    grow.add_argument("--geometry", choices=["surface", "3d", "both"], default="surface", help="Candidate direction mode.")
    grow.add_argument("--directions", type=int, default=8, help="Number of directions sampled per growth centre.")
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
                structure_data, outputs = run_topological_growth(input_file, metals, args.add, args.geometry, args.directions, args.output_dir)
                total_outputs += len(outputs)

                print(
                f"{input_file.name}: "
                f"{structure_data['n_metals']} metal atom(s), "
                f"{len(structure_data['growth_centres'])} centre(s), "
                f"{len(structure_data['growth_edges'])} edge(s), "
                f"{len(outputs)} candidate structure(s)")

            except Exception as error:
                print(f"WARNING: Could not process {input_file}: {error}")

        print(f"Generated {total_outputs} candidate structure(s).")


if __name__ == "__main__":
    main()
