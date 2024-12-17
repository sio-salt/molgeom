import re

from .vec3 import Vec3
from .mat3 import Mat3


def symmop_from_xyz_str(
    xyz_str: str
) -> (Mat3, Vec3):
    """
    Make rotation matrix and translation vector from symmetry operation xyz string.
    Args:
        xyz_str (str): xyz string of symmetry operation.
            e.g.   ‘x, y, z’, ‘-x, -y, z’,
                   '-x, y + 1/2, -z + 1/2', ‘-2y+1/2, 3x+1/2, z-y+1/2’,
    Returns:
        (Mat3, Vec3): Rotation matrix and translation vector in fractional coordinates.
        if you want to convert transvec to cartesian, use molgeom.frac2cart(transvec).
    """

    ops = xyz_str.strip().lower().replace(" ", "").split(",")
    re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
    re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")

    rot_mat = Mat3([[0.0] * 3 for _ in range(3)])
    trans_vec_frac = Vec3(0, 0, 0)
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
