import math

import numpy as np
from numpy.typing import ArrayLike

from .mat3 import Mat3
from .vec3 import Vec3


def lat_params_to_lat_vecs(a, b, c, alpha, beta, gamma, angle_in_degrees=True) -> Mat3:
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

    mat = Mat3(
        [
            a * Vec3(1, 0, 0),
            b * Vec3(cosg, sing, 0),
            c * Vec3(cosb, (cosa - cosb * cosg) / sing, volume / sing),
        ]
    )

    return mat


def lat_vecs_to_lat_params(lattice_vecs: Mat3) -> tuple:
    a = lattice_vecs[0].norm()
    b = lattice_vecs[1].norm()
    c = lattice_vecs[2].norm()
    alpha = math.degrees(lattice_vecs[1].angle(lattice_vecs[2]))
    beta = math.degrees(lattice_vecs[0].angle(lattice_vecs[2]))
    gamma = math.degrees(lattice_vecs[0].angle(lattice_vecs[1]))

    return a, b, c, alpha, beta, gamma


def wrap_frac_coords(frac_coords: np.ndarray) -> np.ndarray:
    return frac_coords % 1.0


def cart2frac(
    cart_coords: ArrayLike, lattice_vecs: ArrayLike, wrap: bool = True
) -> np.ndarray:
    cart_coords = np.asarray(cart_coords)
    lattice_vecs = np.asarray(lattice_vecs)
    if lattice_vecs.shape != (3, 3):
        raise ValueError("lattice_vecs must be a 3x3 matrix")
    if np.linalg.matrix_rank(lattice_vecs) != 3:
        raise ValueError("lattice_vecs must be linearly independent")

    # same as `np.linalg.inv(lattice_vecs).T @ cart_coords` but faster and stable
    frac_coords = np.linalg.solve(lattice_vecs.T, cart_coords)
    if wrap:
        frac_coords = wrap_frac_coords(frac_coords)
    return frac_coords


def frac2cart(frac_coords: Vec3, lattice_vecs: Mat3) -> Vec3:
    return lattice_vecs.T() @ frac_coords
