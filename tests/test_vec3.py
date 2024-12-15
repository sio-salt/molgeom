from molgeom import Vec3
from molgeom import Mat3


def test_vec3_init():
    v1 = Vec3(1, 2, 3)
    assert v1.x == 1 and v1.y == 2 and v1.z == 3


def test_vec3_eq_ne():
    v1 = Vec3(1, 2, 3)
    v2 = Vec3(1, 2, 3)
    v3 = Vec3(4, 5, 6)
    assert v1 == v2 and v1 is not v2
    assert v1 != v3 and v1 is not v3


def test_vec3_add_sub():
    v1 = Vec3(1, 2, 3)
    tmp_id_1 = id(v1)
    v2 = Vec3(4, 5, 6)
    tmp_id_2 = id(v2)
    v3 = v1 + v2
    v4 = v1 - v2
    assert v3 == Vec3(5, 7, 9) and id(v3) != id(v1) and id(v3) != id(v2)
    assert v4 == Vec3(-3, -3, -3) and id(v4) != id(v1) and id(v4) != id(v2)
    v1 += v2
    v2 -= v1
    assert v1 == Vec3(5, 7, 9) and id(v1) == tmp_id_1
    assert v2 == Vec3(-1, -2, -3) and id(v2) == tmp_id_2


def test_vec3_mul_div():
    v1 = Vec3(1, 2, 3)
    v2 = v1 * 2
    v3 = 2 * v1
    v4 = v1 / 2
    v5 = v1 // 2
    assert v2 == Vec3(2, 4, 6) and id(v2) != id(v1)
    assert v3 == Vec3(2, 4, 6) and id(v3) != id(v1)
    assert v4 == Vec3(1 / 2, 2 / 2, 3 / 2) and id(v4) != id(v1)
    assert v5 == Vec3(1 // 2, 2 // 2, 3 // 2) and id(v5) != id(v1)

    v1 = Vec3(1, 2, 3)
    tmp_id_1 = id(v1)
    v1 *= 2
    assert v1 == Vec3(1 * 2, 2 * 2, 3 * 2) and id(v1) == tmp_id_1
    v1 = Vec3(1, 2, 3)
    tmp_id_1 = id(v1)
    v1 /= 2
    assert v1 == Vec3(1 / 2, 2 / 2, 3 / 2) and id(v1) == tmp_id_1
    v1 = Vec3(1, 2, 3)
    tmp_id_1 = id(v1)
    v1 //= 2
    assert v1 == Vec3(1 // 2, 2 // 2, 3 // 2) and id(v1) == tmp_id_1


def test_vec3_neg():
    v1 = Vec3(1, 2, 3)
    v2 = -v1
    assert v2 == Vec3(-1, -2, -3) and id(v2) != id(v1)


def test_vec3_len():
    v1 = Vec3(1, 2, 3)
    assert len(v1) == 3


def test_vec3_dot_cross():
    v1 = Vec3(1, 2, 3)
    v2 = Vec3(4, 5, 6)
    assert v1.dot(v2) == 1 * 4 + 2 * 5 + 3 * 6
    assert v1.cross(v2) == Vec3(2 * 6 - 3 * 5, 3 * 4 - 1 * 6, 1 * 5 - 2 * 4)


def test_vec3_copy():
    v1 = Vec3(1, 2, 3)
    v2 = v1.copy()
    assert v1 == v2 and id(v1) != id(v2)


def test_vec3_str_repr():
    v1 = Vec3(1, 2, 3)
    assert str(v1) == "     1.000000000000      2.000000000000      3.000000000000"
    assert repr(v1) == "Vec3(1.000000000000, 2.000000000000, 3.000000000000)"


def test_vec3_getter_setter():
    v1 = Vec3(1, 2, 3)
    v1.x = 4
    v1.y = 5
    v1.z = 6
    assert v1 == Vec3(4, 5, 6)
    assert v1.x == 4 and v1.y == 5 and v1.z == 6

    try:
        v1.x = "a"
    except TypeError as e:
        assert str(e) == "x must be an int or float"
    try:
        v1.y = "a"
    except TypeError as e:
        assert str(e) == "y must be an int or float"
    try:
        v1.z = "a"
    except TypeError as e:
        assert str(e) == "z must be an int or float"


def test_vec3_types():
    v1 = Vec3(1, 2, 3)
    assert isinstance(v1, Vec3)


def test_mat3_matmul():
    v1 = Vec3(1, 2, 3)
    mat1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    mat2 = Mat3([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    mat3 = [Vec3(1, 2, 3), Vec3(4, 5, 6), Vec3(7, 8, 9)]
    mat4 = Mat3([Vec3(1, 2, 3), Vec3(4, 5, 6), Vec3(7, 8, 9)])

    assert v1.matmul(mat1) == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
    assert v1.matmul(mat2) == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
    assert v1.matmul(mat3) == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
    assert v1.matmul(mat4) == Vec3(
        1 * 1 + 2 * 2 + 3 * 3, 1 * 4 + 2 * 5 + 3 * 6, 1 * 7 + 2 * 8 + 3 * 9
    )
