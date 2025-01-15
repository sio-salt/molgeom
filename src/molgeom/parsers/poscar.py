from __future__ import annotations

import re
from pathlib import Path

from molgeom.utils.vec3 import Vec3
from molgeom.utils.mat3 import Mat3
from molgeom.utils.lattice_utils import frac2cart
from molgeom.atom import Atom
from molgeom.data.consts import ATOMIC_MASSES
from molgeom.molecule import Molecule


def poscar_parser(filepath: str | Path) -> Molecule:
    """
    Parse a POSCAR file and return a Molecule object.

    https://www.vasp.at/wiki/index.php/POSCAR#Full_format_specification

    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"{filepath} do not exist")
    if not filepath.is_file():
        raise ValueError(f"{filepath} is not a file")
    if "poscar" not in filepath.stem.lower():
        raise ValueError(f"{filepath} is not a POSCAR file")

    mol = Molecule()
    with open(filepath, "r") as file:

        # replace tabs, non-breaking spaces, and multiple spaces with single space
        lines = file.readlines()
        for i in range(len(lines)):
            lines[i] = re.sub(r"[\s\t\xa0]+", " ", lines[i])
        lines_gen = iter(lines)

        # The first line is in principle a comment line
        line1 = next(lines_gen).strip()

        # Scaling factor
        scale = list(map(float, next(lines_gen).split()[:3]))
        if len(scale) not in [1, 3]:
            raise RuntimeError("The number of scaling factors must be 1 or 3.")
        if len(scale) == 3 and any([x <= 0 for x in scale]):
            raise RuntimeError("All three scaling factors must be positive.")

        # Lattice vectors
        lattice_vec_a = Vec3(*map(float, next(lines_gen).split()[:3]))
        lattice_vec_b = Vec3(*map(float, next(lines_gen).split()[:3]))
        lattice_vec_c = Vec3(*map(float, next(lines_gen).split()[:3]))
        lattice_vecs = Mat3([lattice_vec_a, lattice_vec_b, lattice_vec_c])
        mol.lattice_vecs = lattice_vecs

        if len(scale) == 1:
            # Negative scaling factor corresponds to the cell volume.
            scale = scale[0]
            if scale < 0.0:
                scale = (-1.0 * scale / lattice_vecs.det()) ** (1 / 3)
            lattice_vecs = Mat3([scale * vec for vec in lattice_vecs])
        else:
            lattice_vecs = Mat3([scale[i] * vec for i, vec in enumerate(lattice_vecs)])

        # Atom symbols and number of atoms per symbol
        atom_symbols = next(lines_gen).split()
        for symbol in atom_symbols:
            if symbol not in ATOMIC_MASSES:
                raise ValueError(f"Unknown atom symbol {symbol=}")
        num_atoms_per_symbol = list(map(int, next(lines_gen).split()))

        # Selective dynamics (optional line)
        tmp_line = next(lines_gen).strip()
        selective_dynamics = False
        if tmp_line.lower()[0] == "s":
            selective_dynamics = True

        # Direct or Cartesian coordinates
        if not selective_dynamics:
            coord_type = tmp_line.strip().lower()
        else:
            coord_type = next(lines_gen).strip().lower()

        # Read the atomic coordinates
        if coord_type == "direct":
            for i in range(len(atom_symbols)):
                for _ in range(num_atoms_per_symbol[i]):
                    frac_coords = Vec3(*map(float, next(lines_gen).split()[:3]))
                    cart_coords = frac2cart(frac_coords, lattice_vecs)
                    atom = Atom(
                        symbol=atom_symbols[i],
                        x=cart_coords[0],
                        y=cart_coords[1],
                        z=cart_coords[2],
                    )
                    mol.add_atom(atom)
        elif coord_type == "cartesian":
            for i in range(len(atom_symbols)):
                for _ in range(num_atoms_per_symbol[i]):
                    cart_coords = list(map(float, next(lines_gen).split()[:3]))
                    atom = Atom(
                        symbol=atom_symbols[i],
                        x=cart_coords[0],
                        y=cart_coords[1],
                        z=cart_coords[2],
                    )
                    mol.add_atom(atom)
        else:
            raise ValueError(f"Expected 'Direct' or 'Cartesian', got {coord_type=}")

    return mol
