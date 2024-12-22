import numpy as np
from numpy.typing import ArrayLike

from .vec3 import Tvec


def lat_params_to_lat_vecs(
    a, b, c, alpha, beta, gamma, angle_in_degrees=True
) -> np.ndarray:
    if angle_in_degrees:
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)

    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = np.sqrt(1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg)

    mat = np.array(
        [
            [a, 0, 0],
            [b * cosg, b * sing, 0],
            [c * cosb, c * ((cosa - cosb * cosg) / sing), c * volume / sing],
        ]
    )

    return mat


def lat_vecs_to_lat_params(
    lattice_vecs: ArrayLike,
) -> tuple[float, float, float, float, float, float]:
    lattice_vecs = np.asarray(lattice_vecs)

    a = np.linalg.norm(lattice_vecs[0])
    b = np.linalg.norm(lattice_vecs[1])
    c = np.linalg.norm(lattice_vecs[2])

    alpha = np.degrees(np.arccos(lattice_vecs[1].dot(lattice_vecs[2]) / (b * c)))
    beta = np.degrees(np.arccos(lattice_vecs[0].dot(lattice_vecs[2]) / (a * c)))
    gamma = np.degrees(np.arccos(lattice_vecs[0].dot(lattice_vecs[1]) / (a * b)))

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

    frac_coords = cart_coords @ np.linalg.inv(lattice_vecs)
    if wrap:
        frac_coords = wrap_frac_coords(frac_coords)
    return frac_coords


def frac2cart(frac_coords: Tvec | ArrayLike, lattice_vecs: ArrayLike) -> np.ndarray:
    frac_coords = np.asarray(frac_coords)
    lattice_vecs = np.asarray(lattice_vecs)
    if lattice_vecs.shape != (3, 3):
        raise ValueError("lattice_vecs must be a 3x3 matrix")
    if np.linalg.matrix_rank(lattice_vecs) != 3:
        raise ValueError("lattice_vecs must be linearly independent")

    return frac_coords @ lattice_vecs
