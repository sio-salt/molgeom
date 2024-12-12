import sys
from collections import deque

from molgeom import Vec3
from molgeom.atom import Atom
from molgeom.molecule import Molecule
from molgeom.parsers.parser_tools import is_valid_xyz_line, remove_trailing_empty_lines


# Gaussian input file parser (com, gjf)
def gau_inp_parser(filepath: str) -> Molecule:
    mole = Molecule()

    with open(filepath, "r") as file:
        # remove trailing empty lines and create a deque to pop from left
        lines = deque(remove_trailing_empty_lines(file.readlines()))

        # link0 section
        link0 = []
        while lines and lines[0].strip().startswith("%"):
            link0.append(lines.popleft().strip())

        # route section
        route_section = []
        if not lines[0].strip().startswith("#"):
            raise ValueError("Expected route section.\n" + f"Found: {lines[0]}\n")
        while lines and lines[0].strip():
            route_section.append(lines.popleft().strip())
        lines.popleft()

        # title section
        if not lines[0].strip() or lines[1].strip():
            raise ValueError(
                "Expected title section.\n"
                + f"Found: {lines[0]}\n"
                + f"Found: {lines[1]}\n"
            )
        _title = lines.popleft().strip()
        lines.popleft()

        # molecule Specification section
        if not lines[0].strip():
            raise ValueError(
                "Expected molecule specification section.\n" + f"Found: {lines[0]}\n"
            )

        try:
            charge, multiplicity = map(int, lines.popleft().strip().split())
        except ValueError as e:
            print(f"gau_parser : invalid file format \n{e}")
            sys.exit(1)

        # atom cartesian coordinates
        while lines and lines[0].strip():
            if not is_valid_xyz_line(lines[0]):
                raise ValueError(
                    "Expected atom symbol and cartesian coordinates.\n"
                    + f"Found: {lines[0]}\n"
                )
            data = lines.popleft().strip().split()
            if data[0] == "Tv":
                if not (
                    lines[0].strip().startswith("Tv")
                    and lines[1].strip().startswith("Tv")
                ):
                    raise ValueError("Expected 3 lines of Tv lattice vectors.")

                lat_vec = [data[1:4]]
                for _ in range(2):
                    if not is_valid_xyz_line(lines[0]):
                        raise ValueError(
                            "Expected atom symbol and cartesian coordinates.\n"
                            + f"Found: {lines[0]}\n"
                        )
                    data = lines.popleft().strip().split()
                    lat_vec.append(data[1:4])
                mole.lattice_vecs = [Vec3(*map(float, vec)) for vec in lat_vec]
            else:
                atom = Atom(
                    symbol=data[0], x=float(data[1]), y=float(data[2]), z=float(data[3])
                )
                mole.add_atom(atom)
        if not mole:
            raise ValueError("No atoms found in file.")

    return mole
