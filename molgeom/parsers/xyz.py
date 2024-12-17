import re

from molgeom.atom import Atom
from molgeom.molecule import Molecule
from molgeom.parsers.parser_tools import (
    atom_xyz_regex,
    is_valid_xyz_line,
    remove_trailing_empty_lines,
)


def xyz_parser(filepath: str) -> Molecule:
    mole = Molecule()

    with open(filepath, "r") as file:
        lines = remove_trailing_empty_lines(file.readlines())

        # replace tabs, non-breaking spaces, and multiple spaces with single space
        for i in range(len(lines)):
            lines[i] = re.sub(r"[\s\t\xa0]+", " ", lines[i])

        first_line = lines[0].strip()
        mol_struc_lines = []
        if re.fullmatch(r"\d+", first_line):
            num_atoms = int(first_line)
            _comment_line = lines[1].strip()
            mol_struc_lines = lines[2:]
        elif re.fullmatch(atom_xyz_regex, first_line):
            mol_struc_lines = lines
            num_atoms = len(mol_struc_lines)
        else:
            raise ValueError(f"Invalid file format in first line: \n{first_line}")

        for line in mol_struc_lines:
            if not is_valid_xyz_line(line):
                raise ValueError(f"Invalid line format: \n{line}")

            data = line.strip().split()
            atom = Atom(
                symbol=data[0], x=float(data[1]), y=float(data[2]), z=float(data[3])
            )
            mole.atoms.append(atom)

        if num_atoms != len(mole):
            raise ValueError(
                f"Number of atoms specified ({num_atoms}) "
                + f"does not match number of atoms read ({len(mole)})."
            )

    return mole
