from __future__ import annotations

import sys
from pathlib import Path

from molgeom.molecule import Molecule
from molgeom.parsers.cif import cif_parser
from molgeom.parsers.gamess import gms_inp_parser
from molgeom.parsers.gaussian import gau_inp_parser
from molgeom.parsers.poscar import poscar_parser
from molgeom.parsers.xyz import xyz_parser
from molgeom.parsers.mol import mol_parser
from molgeom.parsers.sdf import sdf_parser


def read_file(filepath: str | Path) -> Molecule:
    """Read file and return Molecule object."""
    filepath = Path(filepath)

    ext_parser_map = {
        ".xyz": xyz_parser,
        ".com": gau_inp_parser,
        ".gjf": gau_inp_parser,
        ".inp": gms_inp_parser,
        ".cif": cif_parser,
        ".mol": mol_parser,
        ".sdf": sdf_parser,
    }

    if not filepath.exists():
        raise FileNotFoundError(f"{filepath} do not exist")

    for ext in ext_parser_map:
        if filepath.suffix == ext:
            return ext_parser_map[ext](filepath)

    if "poscar" in str(filepath).lower():
        return poscar_parser(filepath)

    raise RuntimeError(
        f'file extension for "{filepath.name}" ' + "is not supported or extensionless file"
    )


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <xyz_file_path>")
        sys.exit(1)

    print()
    filepaths = sys.argv[1:]
    for filepath in filepaths:
        print(filepath)
        mole = read_file(filepath)
        print(mole.to_xyz())


if __name__ == "__main__":
    main()
