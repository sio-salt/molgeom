import math

import numpy as np

from molgeom import Vec3
from molgeom.utils import (
    lat_vecs_to_lat_params,
    lat_params_to_lat_vecs,
    cart2frac,
    frac2cart,
)


def test_lat_vecs_to_lat_params():

    a = 5.0
    b = 5.0
    c = 5.0
    alpha = 90.0
    beta = 90.0
    gamma = 90.0

    lattice_vecs = lat_params_to_lat_vecs(a, b, c, alpha, beta, gamma)

    a2, b2, c2, alpha2, beta2, gamma2 = lat_vecs_to_lat_params(lattice_vecs)

    assert math.isclose(a, a2)
    assert math.isclose(b, b2)
    assert math.isclose(c, c2)
    assert math.isclose(alpha, alpha2)
    assert math.isclose(beta, beta2)
    assert math.isclose(gamma, gamma2)

    params = (1.000000, 2.022375, 7.051241, 84.446654, 85.118719, 81.469234)
    lat_vecs = lat_params_to_lat_vecs(*params)
    tmp_mat = np.asarray(
        [
            Vec3(1.000000, 0.000000, 0.000000),
            Vec3(0.300000, 2.000000, 0.000000),
            Vec3(0.600000, 0.600000, 7.000000),
        ]
    )
    for i in range(3):
        for j in range(3):
            assert math.isclose(lat_vecs[i][j], tmp_mat[i][j], rel_tol=1e-5)

    lat_params_from_tmp = lat_vecs_to_lat_params(tmp_mat)
    for i in range(len(lat_params_from_tmp)):
        assert math.isclose(params[i], lat_params_from_tmp[i], rel_tol=1e-5)


def test_conv_coords():
    lat_vecs = np.asarray(
        [
            [1.000000, 0.000000, 0.000000],
            [0.300000, 2.000000, 0.000000],
            [0.600000, 0.600000, 7.000000],
        ]
    )
    vec1 = Vec3(0.4814285714285711, 0.8714285714285714, 0.4285714285714286)
    vec2 = Vec3(1.0, 2.0, 3.0)
    vec1_cart = frac2cart(vec1, lat_vecs)
    vec2_frac = cart2frac(vec2, lat_vecs)
    assert np.allclose(vec1_cart, vec2)
    assert np.allclose(vec2_frac, vec1)

    cart_coords = [Vec3(0.1, 0.2, 0.3), Vec3(0.4, 0.5, 0.6), Vec3(0.7, 0.8, 0.9)]
    frac_coords = [cart2frac(x, lat_vecs) for x in cart_coords]
    cart_coords2 = [frac2cart(x, lat_vecs) for x in frac_coords]
    for i in range(3):
        for j in range(3):
            assert math.isclose(cart_coords[i][j], cart_coords2[i][j], rel_tol=1e-5)
