import glob
from pathlib import Path
from ase.io import read, write
from .config import default_extensions


def get_input_files(inputs, extensions=None):
    '''
    Expand filenames, directories and glob patterns into a sorted list of input files.
    '''
    if extensions is None:
        extensions = default_extensions

    found = set()

    for raw_input in inputs:
        path = Path(raw_input).expanduser()

        if path.is_file():
            found.add(path.resolve())
            continue

        if path.is_dir():
            for candidate in path.iterdir():
                if candidate.is_file() and candidate.suffix.lower() in extensions:
                    found.add(candidate.resolve())
            continue

        for match in glob.glob(str(path)):
            candidate = Path(match)

            if candidate.is_file():
                found.add(candidate.resolve())
            elif candidate.is_dir():
                for subcandidate in candidate.iterdir():
                    if subcandidate.is_file() and subcandidate.suffix.lower() in extensions:
                        found.add(subcandidate.resolve())

    return sorted(found)


def read_structure(path):
    '''
    Read one structure with ASE.
    '''
    return read(str(path))


def write_structure(path, atoms):
    '''
    Write one ASE structure.
    '''
    write(str(path), atoms)
