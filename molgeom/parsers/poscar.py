from __future__ import annotations
import os
from easyvec import Vec3
from molgeom.data.consts import ATOMIC_MASSES
from molgeom.atom import Atom
from molgeom.molecule import Molecule


def poscar_parser(filepath: str) -> Molecule:
    """
    Parse a POSCAR file and return a Molecule object.

    https://www.vasp.at/wiki/index.php/POSCAR#Full_format_specification

    """
    if not os.path.exists(filepath) or not os.path.isfile(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")
    if "poscar" not in os.path.basename(filepath).lower():
        raise ValueError(f"{filepath} is not a POSCAR file")

    mole = Molecule()
    with open(filepath, "r") as file:
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
        mole.lattice_vecs = [lattice_vec_a, lattice_vec_b, lattice_vec_c]
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

    return mole
