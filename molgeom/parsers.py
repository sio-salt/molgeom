import sys
import os
import re
from collections import deque
from molgeom.data.consts import ATOMIC_MASSES
from molgeom.atom import Atom
from molgeom.molecule import Molecule


atom_xyz_regex = re.compile(r"(\w\w?)(\s+[-+]?\d*\.\d+){3}")


def remove_trailing_empty_lines(lines: list[str]):
    while lines and not lines[-1].strip():
        lines.pop()
    return lines


def is_valid_xyz_line(line: str) -> bool:
    data = line.strip().split()
    if not data:
        return False
    if len(data) != 4:
        return False
    if not re.fullmatch(atom_xyz_regex, line.strip()):
        return False
    if data[0] not in ATOMIC_MASSES:
        return False
    return True


def xyz_parser(filepath: str, mode: str) -> Molecule:
    mole = Molecule()

    with open(filepath, mode) as file:
        lines = remove_trailing_empty_lines(file.readlines())

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
            print("xyz_parser : invalid file format in first line \n" + f"{first_line}")
            sys.exit(1)

        for line in mol_struc_lines:
            if not is_valid_xyz_line(line):
                print(
                    "xyz_parser : invalid data format or empty line found " + f"{line}"
                )
                sys.exit(1)

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


# Gaussian com file parser
def com_parser(filepath: str, mode: str) -> Molecule:
    mole = Molecule()

    with open(filepath, mode) as file:
        lines = deque(remove_trailing_empty_lines(file.readlines()))

        # link0 section
        link0 = []
        while lines and lines[0].strip().startswith("%"):
            link0.append(lines.popleft().strip())

        # route section
        route_section = []
        if not lines[0].strip().startswith("#"):
            print(f"com_parser : invalid file format \n{lines[0]}")
            sys.exit(1)
        while lines and lines[0].strip():
            route_section.append(lines.popleft().strip())
        lines.popleft()

        # title section
        if not lines[0].strip() or lines[1].strip():
            print(f"com_parser : invalid file format \n{lines[0]}")
            sys.exit(1)
        _title = lines.popleft().strip()
        lines.popleft()

        # molecule Specification section
        if not lines[0].strip():
            print(f"com_parser : invalid file format \n{lines[0]}")
            sys.exit(1)
        try:
            charge, multiplicity = map(int, lines.popleft().strip().split())
        except ValueError as e:
            print(f"com_parser : invalid file format \n{e}")
            sys.exit(1)

        # atom cartesian coordinates
        while lines and lines[0].strip():
            if not is_valid_xyz_line(lines[0]):
                print(f"com_parser : invalid file format \n{lines[0]}")
                sys.exit(1)
            data = lines.popleft().strip().split()
            atom = Atom(
                symbol=data[0], x=float(data[1]), y=float(data[2]), z=float(data[3])
            )
            mole.add_atoms(atom)
        if not mole:
            raise ValueError("No atoms found in file.")

    return mole


atom_mass_xyz_regex = re.compile(r"(\w\w?)(\s+\d+\.\d+)(\s+[-+]?\d*\.\d+){3}")


def is_valid_gms_xyz_line(line: str) -> bool:
    data = line.strip().split()
    if not data:
        return False
    if len(data) != 5:
        return False
    if not re.fullmatch(atom_mass_xyz_regex, line.strip()):
        return False
    if data[0] not in ATOMIC_MASSES:
        return False
    return True


# GAMESS input file parser
def inp_parser(filepath: str, mode: str) -> Molecule:
    mole = Molecule()

    with open(filepath, mode) as file:
        lines = deque(remove_trailing_empty_lines(file.readlines()))

        # input description section
        input_descriptions = []
        while lines and not lines[0].strip().upper().startswith("$DATA"):
            if lines[0].strip().startswith("!"):
                lines.popleft()
                continue
            input_descriptions.append(lines.popleft().strip())

        # $DATA group
        lines.popleft()
        if not lines[0].strip() or not lines[1].strip():
            print(
                "inp_parser DATA group : invalid file format \n"
                + f"{lines[0]}\n"
                + f"{lines[1]}"
            )
            sys.exit(1)
        _title = lines.popleft().strip()
        _group_naxis = lines.popleft().strip()

        # atom cartesian coordinates
        while lines and lines[0].strip() != "$END":
            if not is_valid_gms_xyz_line(lines[0]):
                print(
                    "inp_parser atom cartesian coords : "
                    + f"invalid file format \n{lines[0]}"
                )
                sys.exit(1)
            data = lines.popleft().strip().split()
            atom = Atom(
                symbol=data[0], x=float(data[2]), y=float(data[3]), z=float(data[4])
            )
            mole.add_atoms(atom)

    return mole


def parse_file(filepath: str, mode: str = "r") -> Molecule:
    ext_parser_map = {".xyz": xyz_parser, ".com": com_parser, ".inp": inp_parser}
    supported_mode = {"r"}

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")
    if mode not in supported_mode:
        raise ValueError(f"mode {mode} is not supported yet")

    for ext in ext_parser_map:
        if filepath.endswith(ext):
            return ext_parser_map[ext](filepath, mode)
    else:
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
        mole = parse_file(filepath, "r")
        print(mole.to_xyz())


if __name__ == "__main__":
    main()
