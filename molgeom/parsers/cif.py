import os
import math
from collections import deque

from molgeom import Vec3, Atom, Molecule
from molgeom.parsers.parser_tools import remove_trailing_empty_lines

# CIF format specification:
# http://www.physics.gov.az/book_I/S_R_Hall.pdf


def cif_tag_parser(filepath: str) -> dict:
    if not os.path.exists(filepath) or not os.path.isfile(filepath):
        raise FileNotFoundError(f"{filepath} do not exist")

    cif_tags = dict()
    with open(filepath, "r") as file:
        lines = deque()
        # add empty line before loop_ block to separate tags
        for line in remove_trailing_empty_lines(file.readlines()):
            if "loop_" in line:
                lines.append(" ")
                lines.append(line)
            else:
                lines.append(line)
        lines.append(" ")
        lines.append(" ")

        while lines:
            line = lines.popleft()

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

            # get tags under loop_ block
            loop_tags_idx = dict()
            if line.strip().startswith("loop_"):
                line = lines.popleft().strip()
                cnt = 0
                while line.startswith("_"):
                    tag = line.split()[0]
                    loop_tags_idx[tag] = cnt
                    line = lines.popleft().strip()
                    cnt += 1

            # load symmetry operations
            symop_tags = [
                "_space_group_symop_operation_xyz",
                "_space_group_symop.operation_xyz",
                "_symmetry_equiv_pos_as_xyz",
            ]
            symop_tags_used = [tag for tag in loop_tags_idx.keys() if tag in symop_tags]
            if symop_tags_used:
                symop_idx = loop_tags_idx[symop_tags_used[0]]
                cif_tags["symops"] = []
                while line.count("'") == 2:
                    quote_idxes = [
                        line.index("'"),
                        line.index("'", line.index("'") + 1),
                    ]
                    symop_str = line[quote_idxes[0] + 1 : quote_idxes[1]]
                    if symop_str.count(",") != 2:
                        raise ValueError(
                            f"Invalid symop string: {symop_str}\n"
                            + "at least 3 tokens separated by comma are required"
                        )
                    cif_tags["symops"].append(symop_str)
                    line = lines.popleft().strip()

            # load atom symbols and positions
            atom_tags = [
                "_atom_site_type_symbol",
                "_atom_site_fract_x",
                "_atom_site_fract_y",
                "_atom_site_fract_z",
            ]
            if any(atom_tag in loop_tags_idx for atom_tag in atom_tags):
                if not all(atom_tag in loop_tags_idx for atom_tag in atom_tags):
                    raise ValueError("All atom tags must be present in the loop_ block")
                cif_tags["atoms"] = []
                while len(line.split()) == len(loop_tags_idx):
                    splited = line.split()
                    symbol = splited[loop_tags_idx["_atom_site_type_symbol"]]
                    fract_x = splited[loop_tags_idx["_atom_site_fract_x"]].split("(")[0]
                    fract_y = splited[loop_tags_idx["_atom_site_fract_y"]].split("(")[0]
                    fract_z = splited[loop_tags_idx["_atom_site_fract_z"]].split("(")[0]
                    atom = {
                        "symbol": symbol,
                        "fract_x": float(fract_x),
                        "fract_y": float(fract_y),
                        "fract_z": float(fract_z),
                    }
                    cif_tags["atoms"].append(atom)
                    line = lines.popleft().strip()

    if "atoms" not in cif_tags or len(cif_tags["atoms"]) == 0:
        raise ValueError("No atoms found in the CIF file")
    if any(
        tag not in cif_tags
        for tag in [
            "cell_length_a",
            "cell_length_b",
            "cell_length_c",
            "cell_angle_alpha",
            "cell_angle_beta",
            "cell_angle_gamma",
        ]
    ):
        raise ValueError("Cell parameters not found in the CIF file")
    if "symops" in cif_tags and len(cif_tags["symops"]) == 0:
        raise ValueError("No symmetry operations found in the CIF file")

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
