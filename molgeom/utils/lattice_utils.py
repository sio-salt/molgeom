import math
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


def wrap_frac_coords(frac_coords: Vec3) -> Vec3:
    return frac_coords % 1.0


def cart2frac(cart_coords: Vec3, lattice_vecs: Mat3, wrap=False) -> Vec3:
    frac_coords = lattice_vecs.inv().T() @ cart_coords
    if wrap:
        frac_coords = wrap_frac_coords(frac_coords)
    return frac_coords


def frac2cart(frac_coords: Vec3, lattice_vecs: Mat3) -> Vec3:
    return lattice_vecs.T() @ frac_coords
