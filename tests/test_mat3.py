import math

from molgeom import Vec3
from molgeom import Mat3


def test_mat3_init():
    mat1 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert (
        mat1[0] == Vec3(1, 2, 3)
        and mat1[1] == Vec3(4, 5, 6)
        and mat1[2] == Vec3(7, 8, 9)
    )

    mat1[0] = Vec3(2, 0, -1)
    assert mat1[0] == Vec3(2, 0, -1)
    mat1[0] = [1, 2, 3]
    assert mat1[0] == Vec3(1, 2, 3)
    try:
        mat1[0] = [1, 2, 3, 4]
    except ValueError as e:
        assert str(e) == "Invalid vector type"

    try:
        mat2 = Mat3([1, 2, 3], [4, 5, 6], [7, 8])
    except ValueError as e:
        assert str(e) == "Matrix must be 3x3"
    try:
        mat3 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9, 10])
    except ValueError as e:
        assert str(e) == "Matrix must be 3x3"
    try:
        mat4 = Mat3([1, 2, 3], [4, 5, 6])
    except ValueError as e:
        assert str(e) == "Matrix must be 3x3"


def test_mat3_str_repr():
    mat1 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    assert str(mat1) == "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]"
    assert repr(mat1) == "Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])"


def test_mat3_eq():
    mat1 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    mat2 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    mat3 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 10])
    assert mat1 == mat2 and mat1 != mat3


def test_mat3_add_sub():
    mat1 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    mat2 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    mat3 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 10])

    assert mat1 + mat2 == Mat3([2, 4, 6], [8, 10, 12], [14, 16, 18])
    assert mat1 - mat2 == Mat3([0, 0, 0], [0, 0, 0], [0, 0, 0])
    assert mat1 + mat3 == Mat3([2, 4, 6], [8, 10, 12], [14, 16, 19])
    assert mat1 - mat3 == Mat3([0, 0, 0], [0, 0, 0], [0, 0, -1])

    mat1 += mat2
    assert mat1 == Mat3([2, 4, 6], [8, 10, 12], [14, 16, 18])
    mat1 -= mat2
    assert mat1 == Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])

    try:
        mat1 + 1
    except ValueError as e:
        assert str(e) == "Invalid matrix type"
    try:
        mat1 - 1
    except ValueError as e:
        assert str(e) == "Invalid matrix type"


def test_mat3_mul_div():
    mat1 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])

    assert mat1 * 2 == Mat3([2, 4, 6], [8, 10, 12], [14, 16, 18])
    assert mat1 / 2 == Mat3([0.5, 1, 1.5], [2, 2.5, 3], [3.5, 4, 4.5])

    mat1 *= 2
    assert mat1 == Mat3([2, 4, 6], [8, 10, 12], [14, 16, 18])
    mat1 /= 2
    assert mat1 == Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])

    try:
        mat1 * "a"
    except ValueError as e:
        assert str(e) == "Invalid matrix type"
    try:
        mat1 / "a"
    except ValueError as e:
        assert str(e) == "Invalid matrix type"


def test_mat3_pow():
    mat1 = Mat3([1, 2, 3], [4, 5, 6], [7, 8, 9])
    assert mat1**2 == Mat3([1, 4, 9], [16, 25, 36], [49, 64, 81])
    assert mat1**0 == Mat3([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert mat1**1 == mat1

    mat1 **= 2
    assert mat1 == Mat3([[1, 4, 9], [16, 25, 36], [49, 64, 81]])

    try:
        mat1 ** "a"
    except ValueError as e:
        assert str(e) == "Invalid matrix type"


def test_mat3_iter():
    mat1 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    for i, vec in enumerate(mat1):
        assert vec == mat1[i]
        assert isinstance(vec, Vec3)


def test_mat3_types():
    mat1 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert isinstance(mat1, Mat3)


def test_mat3_matmul():
    v1 = Vec3(1, 2, 3)
    mat1 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    mat2 = Mat3([Vec3(1, 2, 3), Vec3(4, 5, 6), Vec3(7, 8, 9)])
    mat3 = Mat3(
        [[math.cos(1), -math.sin(1), 0], [math.sin(1), math.cos(1), 0], [0, 0, 1]]
    )
    mat4 = Mat3(
        [[math.cos(2), -math.sin(2), 0], [math.sin(2), math.cos(2), 0], [0, 0, 1]]
    )
    mat5 = Mat3(
        [[math.cos(3), -math.sin(3), 0], [math.sin(3), math.cos(3), 0], [0, 0, 1]]
    )

    assert mat1 @ v1 == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
    assert mat2 @ v1 == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
    assert mat1 @ mat2 == Mat3([[30, 36, 42], [66, 81, 96], [102, 126, 150]])
    mat6 = mat3 @ mat4
    for i in range(3):
        for j in range(3):
            assert math.isclose(mat6[i][j], mat5[i][j], rel_tol=1e-8, abs_tol=1e-8)


def test_mat3_matrix_power():
    mat1 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert mat1.matrix_power(2) == mat1 @ mat1
    assert mat1.matrix_power(0) == Mat3([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert mat1.matrix_power(1) == mat1
    assert mat1.matrix_power(3) == mat1 @ mat1 @ mat1


def test_mat3_mat_det():
    identity = Mat3([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert identity.det() == 1

    mat1 = Mat3([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    # det A = A11A22A33 + A12A23A31 + A13A21A32 - A13A22A31 - A12A21A33 - A11A23A32
    mat1_det = (
        (2 * 2 * 2)
        + (-1 * -1 * 0)
        + (0 * -1 * -1)
        - (0 * 2 * 0)
        - (-1 * -1 * 2)
        - (2 * -1 * -1)
    )
    assert mat1.det() == mat1_det

    mat2 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    mat2_det = (
        (1 * 5 * 9)
        + (2 * 6 * 7)
        + (3 * 4 * 8)
        - (3 * 5 * 7)
        - (2 * 4 * 9)
        - (1 * 6 * 8)
    )
    assert mat2.det() == mat2_det

    mat3 = Mat3([[2, 2, 1], [1, -8, 0], [3, -7, 1]])
    mat3_det = (
        (2 * -8 * 1)
        + (2 * 0 * 3)
        + (1 * 1 * -7)
        - (1 * -8 * 3)
        - (2 * 1 * 1)
        - (2 * 0 * -7)
    )
    assert mat3.det() == mat3_det


def test_mat3_inv_mat():
    identity = Mat3([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert identity.inv() == identity

    mat1 = Mat3([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    assert mat1.inv() == [
        [3 / 4, 2 / 4, 1 / 4],
        [2 / 4, 4 / 4, 2 / 4],
        [1 / 4, 2 / 4, 3 / 4],
    ]

    # use sin cos
    mat2 = Mat3(
        [[math.cos(1), -math.sin(1), 0], [math.sin(1), math.cos(1), 0], [0, 0, 1]]
    )
    assert mat2.inv() == Mat3(
        [
            [math.cos(1), math.sin(1), 0],
            [-math.sin(1), math.cos(1), 0],
            [0, 0, 1],
        ]
    )

    mat3 = Mat3([[2, 2, 1], [1, -8, 0], [3, -7, 1]])
    mat3_inv = Mat3(
        [
            [8, 9, -8],
            [1, 1, -1],
            [-17, -20, 18],
        ]
    )
    assert mat3.inv() == mat3_inv
    assert mat3_inv.inv() == mat3
    assert mat3.det() == 1
    assert mat3.det() == 1
