import math
import os

from molgeom import Vec3, Atom, Molecule

# CIF format specification:
# http://www.physics.gov.az/book_I/S_R_Hall.pdf


def cif_tag_parser(filepath: str) -> dict:
    if not os.path.exists(filepath) or not os.path.isfile(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")

    cif_tags = dict()
    with open(filepath, "r") as file:
        for line in file:
            if any(
                line.strip().startswith(tag)
                for tag in [
                    "_space_group_symop_operation_xyz",
                    "_space_group_symop.operation_xyz",
                    "_symmetry_equiv_pos_as_xyz",
                ]
            ):
                line = next(file)
                symops = []
                while line.strip():
                    symop_str = line.strip().replace("'", "").replace('"', "")
                    if symop_str.count(",") != 2:
                        raise ValueError(
                            f"Invalid symop string: {symop_str}\n"
                            + "at least 3 tokens separated by comma are required"
                        )
                    symops.append(symop_str)
                    line = next(file)
                cif_tags["symops"] = symops

            # load cell parameters
            if line.strip().startswith("_cell_length_a"):
                len_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_length_a"] = float(len_val_str)
            if line.strip().startswith("_cell_length_b"):
                len_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_length_b"] = float(len_val_str)
            if line.strip().startswith("_cell_length_c"):
                len_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_length_c"] = float(len_val_str)
            if line.strip().startswith("_cell_angle_alpha"):
                angle_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_angle_alpha"] = float(angle_val_str)
            if line.strip().startswith("_cell_angle_beta"):
                angle_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_angle_beta"] = float(angle_val_str)
            if line.strip().startswith("_cell_angle_gamma"):
                angle_val_str = line.split()[1].split("(")[0]
                cif_tags["cell_angle_gamma"] = float(angle_val_str)

            # load atom symbols and positions
            if line.strip().startswith("_atom_site_label"):
                while line.strip().startswith("_atom"):
                    line = next(file)
                atoms = []
                while line.strip():
                    line_splited = line.strip().split()
                    symbol = line_splited[1]
                    fract_x = line_splited[2].split("(")[0]
                    fract_y = line_splited[3].split("(")[0]
                    fract_z = line_splited[4].split("(")[0]
                    atom = {
                        "symbol": symbol,
                        "fract_x": float(fract_x),
                        "fract_y": float(fract_y),
                        "fract_z": float(fract_z),
                    }
                    atoms.append(atom)
                    line = next(file)
                cif_tags["atoms"] = atoms

    return cif_tags


def get_fract_to_cart_mat(
    len_a, len_b, len_c, alpha, beta, gamma, angle_in_degrees=True
) -> list[list[float]]:
    # rewrite this function and don't use numpy
    if angle_in_degrees:
        alpha = math.radians(alpha)
        beta = math.radians(beta)
        gamma = math.radians(gamma)

    cosa = math.cos(alpha)
    cosb = math.cos(beta)
    cosg = math.cos(gamma)
    sing = math.sin(gamma)
    volume = math.sqrt(
        1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    )

    mat = [
        [len_a, len_b * cosg, len_c * cosb],
        [0, len_b * sing, len_c * (cosa - cosb * cosg) / sing],
        [0, 0, len_c * volume / sing],
    ]

    return mat


def ciftag2mol(cif_tags: dict) -> Molecule:
    frac_to_cart_mat = get_fract_to_cart_mat(
        cif_tags["cell_length_a"],
        cif_tags["cell_length_b"],
        cif_tags["cell_length_c"],
        cif_tags["cell_angle_alpha"],
        cif_tags["cell_angle_beta"],
        cif_tags["cell_angle_gamma"],
    )
    mol = Molecule()
    for atom in cif_tags["atoms"]:
        fract_vec = Vec3(atom["fract_x"], atom["fract_y"], atom["fract_z"])
        cart_vec = fract_vec.matmul(frac_to_cart_mat)
        symbol = atom["symbol"]
        atom = Atom.from_vec(symbol, cart_vec)
        mol.add_atom(atom)

    lattice_vecs = []
    for i in range(3):
        vec = []
        for j in range(3):
            vec.append(frac_to_cart_mat[j][i])
        lattice_vecs.append(Vec3(*vec))
    mol.lattice_vecs = lattice_vecs

    return mol


def cif_parser(filepath: str, apply_symop: bool = True) -> Molecule:
    cif_tags = cif_tag_parser(filepath)
    mol = ciftag2mol(cif_tags)
    tmp_mol = mol.copy()
    if apply_symop and "symops" in cif_tags:
        for symop in cif_tags["symops"]:
            new_mol = tmp_mol.replicated_from_xyz_str(symop)
            mol.merge(new_mol)
    mol.lattice_vecs = tmp_mol.lattice_vecs
    return mol
