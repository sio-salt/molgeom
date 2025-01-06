import os
import re
from collections import deque

from molgeom import Vec3, Atom, Molecule, Mat3
from molgeom.parsers.parser_tools import remove_trailing_empty_lines
from molgeom.utils.lattice_utils import lat_params_to_lat_vecs

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
            # replace tabs, non-breaking spaces, and multiple spaces with single space
            line = re.sub(r"[\s\t\xa0]+", " ", line)
            if "loop_" in line:
                lines.append(" ")
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
                while line.count(",") == 2:
                    line_str = "".join(line.split()[symop_idx:])
                    symop_str = line_str.replace("'", "")
                    cif_tags["symops"].append(symop_str)
                    line = lines.popleft().strip()

            # load atom symbols and positions
            atom_symbol_tags = [
                "_atom_site_type_symbol",
                "_atom_site_label",
            ]
            atom_fract_tags = [
                "_atom_site_fract_x",
                "_atom_site_fract_y",
                "_atom_site_fract_z",
            ]
            if (
                len(
                    {
                        atom_symbol_tag
                        for atom_symbol_tag in atom_symbol_tags
                        if atom_symbol_tag in loop_tags_idx
                    }
                )
                > 0
            ) and (
                len(
                    {
                        atom_tag
                        for atom_tag in atom_fract_tags
                        if atom_tag in loop_tags_idx
                    }
                )
                == 3
            ):
                cif_tags["atoms"] = []
                while len(line.split()) == len(loop_tags_idx):
                    splited = line.split()
                    if atom_symbol_tags[0] in loop_tags_idx:
                        symbol = splited[loop_tags_idx["_atom_site_type_symbol"]]
                        # remove charge info from symbol
                        symbol = re.sub(r"\d+[+-]?", "", symbol)
                    else:
                        symbol = splited[loop_tags_idx["_atom_site_label"]]
                        # remove label info from symbol
                        symbol = re.sub(r"\d+[+-]?", "", re.split("_", symbol)[0])
                        symbol = symbol.replace("HW", "H").replace("OW", "O")
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


def ciftag2mol(cif_tags: dict) -> Molecule:
    frac_to_cart_mat: Mat3 = lat_params_to_lat_vecs(
        cif_tags["cell_length_a"],
        cif_tags["cell_length_b"],
        cif_tags["cell_length_c"],
        cif_tags["cell_angle_alpha"],
        cif_tags["cell_angle_beta"],
        cif_tags["cell_angle_gamma"],
        angle_in_degrees=True,
    )
    mol = Molecule()
    for atom in cif_tags["atoms"]:
        fract_vec = Vec3(atom["fract_x"], atom["fract_y"], atom["fract_z"])
        cart_vec = frac_to_cart_mat @ fract_vec
        symbol = atom["symbol"]
        atom = Atom.from_vec(symbol, cart_vec)
        mol.add_atom(atom)

    mol.lattice_vecs = frac_to_cart_mat

    return mol


def cif_parser(filepath: str, apply_symop: bool = True) -> Molecule:
    cif_tags = cif_tag_parser(filepath)
    mol = ciftag2mol(cif_tags)
    rep_mol = Molecule()
    if apply_symop and "symops" in cif_tags:
        for symop in cif_tags["symops"]:
            new_mol = mol.replicated_from_xyz_str(symop, wrap=True)
            rep_mol.merge(new_mol)
    rep_mol.lattice_vecs = mol.lattice_vecs
    return rep_mol
