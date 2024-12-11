import sys
import os
from pathlib import Path
from molgeom.molecule import Molecule
from molgeom.parsers.xyz import xyz_parser
from molgeom.parsers.gaussian import gau_inp_parser
from molgeom.parsers.gamess import gms_inp_parser
from molgeom.parsers.cif import cif_parser
from molgeom.parsers.poscar import poscar_parser


def read_file(filepath: str | Path) -> Molecule:
    ext_parser_map = {
        ".xyz": xyz_parser,
        ".com": gau_inp_parser,
        ".gjf": gau_inp_parser,
        ".inp": gms_inp_parser,
        ".cif": cif_parser,
    }

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")

    for ext in ext_parser_map:
        if str(filepath).endswith(ext):
            return ext_parser_map[ext](filepath)

    if "poscar" in str(filepath).lower():
        return poscar_parser(filepath)

    raise RuntimeError(
        f'file extension for "{os.path.basename(filepath)}" '
        + "is not supported or extensionless file"
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
