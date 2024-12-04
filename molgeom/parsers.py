import sys
import os
import re
from pathlib import Path
from collections import deque
from easyvec import Vec3
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
            mole.add_atoms_from(atom)
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
            mole.add_atoms_from(atom)

    return mole


def poscar_parser(
    filepath: str,
    mode: str = "r",
    rep_cell_a: list[int] = [0, 1],
    rep_cell_b: list[int] = [0, 1],
    rep_cell_c: list[int] = [0, 1],
) -> Molecule:
    """
    Parse a POSCAR file and return a Molecule object.

    https://www.vasp.at/wiki/index.php/POSCAR#Full_format_specification

    """
    if mode != "r":
        raise ValueError("mode should be 'r'")
    if not os.path.exists(filepath) or not os.path.isfile(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")
    if "poscar" not in os.path.basename(filepath).lower():
        raise ValueError(f"{filepath} is not a POSCAR file")
    mole = Molecule()
    with open(filepath, mode) as file:
        # The first line is in principle a comment line
        line1 = file.readline().strip()

        # Scaling factor
        scale = list(map(float, file.readline().split()[:3]))
        if len(scale) not in [1, 3]:
            raise RuntimeError("The number of scaling factors must be 1 or 3.")
        if len(scale) == 3 and any([x <= 0 for x in scale]):
            raise RuntimeError("All three scaling factors must be positive.")

        # Lattice vectors
        lattice_vec_a = Vec3(*map(float, file.readline().split()[:3]))
        lattice_vec_b = Vec3(*map(float, file.readline().split()[:3]))
        lattice_vec_c = Vec3(*map(float, file.readline().split()[:3]))
        lattice_vecs = [lattice_vec_a, lattice_vec_b, lattice_vec_c]

        if len(scale) == 1:
            # Negative scaling factor corresponds to the cell volume.
            scale = scale[0]
            if scale < 0.0:
                lattice_vecs_det = sum(
                    lattice_vecs[0][i]
                    * (
                        lattice_vecs[1][(i + 1) % 3] * lattice_vecs[2][(i + 2) % 3]
                        - lattice_vecs[1][(i + 2) % 3] * lattice_vecs[2][(i + 1) % 3]
                    )
                    for i in range(3)
                )
                scale = (-1.0 * scale / lattice_vecs_det) ** (1 / 3)
            lattice_vecs = [scale * vec for vec in lattice_vecs]
        else:
            lattice_vecs = [scale[i] * vec for i, vec in enumerate(lattice_vecs)]

        # Atom symbols and number of atoms per symbol
        atom_symbols = file.readline().split()
        for symbol in atom_symbols:
            if symbol not in ATOMIC_MASSES:
                raise ValueError(f"Unknown atom symbol {symbol=}")
        num_atoms_per_symbol = list(map(int, file.readline().split()))

        # Selective dynamics (optional line)
        tmp_line = file.readline().strip()
        selective_dynamics = False
        if tmp_line.lower()[0] == "s":
            selective_dynamics = True

        # Direct or Cartesian coordinates
        if not selective_dynamics:
            coord_type = tmp_line.strip().lower()
        else:
            coord_type = file.readline().strip().lower()

        # Read the atomic coordinates
        if coord_type == "direct":
            for i in range(len(atom_symbols)):
                for _ in range(num_atoms_per_symbol[i]):
                    frac_coords = list(map(float, file.readline().split()[:3]))
                    cart_coords = [
                        sum([frac_coords[j] * lattice_vecs[j][i] for j in range(3)])
                        for i in range(3)
                    ]
                    atom = Atom(
                        symbol=atom_symbols[i],
                        x=cart_coords[0],
                        y=cart_coords[1],
                        z=cart_coords[2],
                    )
                    mole.add_atom(atom)
        elif coord_type == "cartesian":
            for i in range(len(atom_symbols)):
                for _ in range(num_atoms_per_symbol[i]):
                    cart_coords = list(map(float, file.readline().split()[:3]))
                    atom = Atom(
                        symbol=atom_symbols[i],
                        x=cart_coords[0],
                        y=cart_coords[1],
                        z=cart_coords[2],
                    )
                    mole.add_atom(atom)
        else:
            raise ValueError(f"Unknown coordinate type {coord_type=}")

    # Replicate the cell
    if (
        rep_cell_a[0] >= rep_cell_a[1]
        or rep_cell_b[0] >= rep_cell_b[1]
        or rep_cell_c[0] >= rep_cell_c[1]
    ):
        raise ValueError(
            "The lower bound of the range must be less than the upper bound."
        )

    result_mole = Molecule()
    for i in range(rep_cell_a[0], rep_cell_a[1]):
        for j in range(rep_cell_b[0], rep_cell_b[1]):
            for k in range(rep_cell_c[0], rep_cell_c[1]):
                tmp_mole = mole.copy()
                trans_vec = Vec3(
                    i * lattice_vecs[0][0]
                    + j * lattice_vecs[1][0]
                    + k * lattice_vecs[2][0],
                    i * lattice_vecs[0][1]
                    + j * lattice_vecs[1][1]
                    + k * lattice_vecs[2][1],
                    i * lattice_vecs[0][2]
                    + j * lattice_vecs[1][2]
                    + k * lattice_vecs[2][2],
                )
                tmp_mole.translate(trans_vec)
                result_mole.merge(tmp_mole)

    return result_mole


def parse_file(filepath: str | Path, mode: str = "r") -> Molecule:
    ext_parser_map = {".xyz": xyz_parser, ".com": com_parser, ".inp": inp_parser}
    supported_mode = {"r"}

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")
    if mode not in supported_mode:
        raise ValueError(f"mode {mode} is not supported yet")

    for ext in ext_parser_map:
        if str(filepath).endswith(ext):
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
