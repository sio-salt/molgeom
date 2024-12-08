#!/usr/bin/python3
import sys
import os
from molgeom import poscar_parser


help_msg = """
Usage:
python poscar2xyz.py <poscar_file_path> a_ini a_fin b_ini b_fin c_ini c_fin

a_ini, a_fin, b_ini, b_fin, c_ini, c_fin: range of cell repetitions. accepts only integers.
e.g. `python poscar2xyz.py POSCAR -1 2 -1 2 -1 2` will generate 3x3x3 supercell with the original cell at the center.
"""
    # rep_cell_a: list[int] = [0, 1],
    # rep_cell_b: list[int] = [0, 1],
    # rep_cell_c: list[int] = [0, 1],

def main():
    if len(sys.argv) < 8:
        print(help_msg)
        sys.exit(1)

    filepath, *cell_rep = sys.argv[1:]
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        print(help_msg)
        sys.exit(1)

    try:
        cell_rep = list(map(int, cell_rep))
    except ValueError:
        print("Invalid cell repetition range. Must be integers.")
        print(help_msg)
        sys.exit(1)

    mole = poscar_parser(filepath)
    mole.replicate(
        [cell_rep[0], cell_rep[1]],
        [cell_rep[2], cell_rep[3]],
        [cell_rep[4], cell_rep[5]]
    )
    print(mole.to_xyz())


if __name__ == "__main__":
    main()
