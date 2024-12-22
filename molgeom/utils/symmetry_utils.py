import re

import numpy as np


def symmop_from_xyz_str(xyz_str: str) -> (np.ndarray, np.ndarray):
    """
    Make rotation matrix and translation vector from symmetry operation xyz string.
    Args:
        xyz_str (str): xyz string of symmetry operation.
            e.g.   ‘x, y, z’, ‘-x, -y, z’,
                   '-x, y + 1/2, -z + 1/2', ‘-2y+1/2, 3x+1/2, z-y+1/2’,
    Returns:
        (np.ndarray, np.ndarray): (3, 3) and (3,) shaped numpy arrays in fractional coordinates.
        if you want to convert transvec to cartesian, use molgeom.frac2cart(transvec).
    """

    ops = xyz_str.strip().lower().replace(" ", "").split(",")
    re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
    re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")

    rot_mat = np.zeros((3, 3))
    trans_vec_frac = np.zeros(3)
    for i, op in enumerate(ops):

        # make rot mat
        for match in re_rot.finditer(op):
            # match[0] contains the whole match
            # match[n] (n > 0) contains the n-th group surrounded by ()
            factor = -1.0 if match[1] == "-" else 1.0
            if match[2] != "":
                factor *= (
                    float(match[2]) / float(match[3])
                    if match[3] != ""
                    else float(match[2])
                )
            j = ord(match[4]) - 120
            rot_mat[i][j] = factor

        # make trans vec
        for match in re_trans.finditer(op):
            factor = -1 if match[1] == "-" else 1
            num = (
                float(match[2]) / float(match[3]) if match[3] != "" else float(match[2])
            )
            trans_vec_frac[i] = factor * num

    return rot_mat, trans_vec_frac
