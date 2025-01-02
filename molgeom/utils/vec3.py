from __future__ import annotations
from copy import deepcopy

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .mat3 import Mat3

# import numpy as np
# from scipy.spatial.transform import Rotation as R


class Vec3:
    def __init__(self, x: int | float, y: int | float, z: int | float) -> None:
        if not all(isinstance(i, (int, float)) for i in (x, y, z)):
            raise TypeError("x, y, z must be an int or float")
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    @property
    def x(self) -> float:
        return self._x

    @x.setter
    def x(self, value: int | float) -> None:
        if not isinstance(value, (int, float)):
            raise TypeError("x must be an int or float")
        self._x = float(value)

    @property
    def y(self) -> float:
        return self._y

    @y.setter
    def y(self, value: int | float) -> None:
        if not isinstance(value, (int, float)):
            raise TypeError("y must be an int or float")
        self._y = float(value)

    @property
    def z(self) -> float:
        return self._z

    @z.setter
    def z(self, value: int | float) -> None:
        if not isinstance(value, (int, float)):
            raise TypeError("z must be an int or float")
        self._z = float(value)

    @classmethod
    def from_list(cls, vec: list[int | float]) -> Vec3:
        if not is_vec_type(vec):
            raise ValueError("vec must be a list of 3 int or float")
        return cls(*vec)

    def to_list(self) -> list[float]:
        return [self.x, self.y, self.z]

    def to_dict(self) -> dict[str, float]:
        return {"x": self.x, "y": self.y, "z": self.z}

    def __str__(self) -> str:
        return f"{self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    def __repr__(self) -> str:
        return f"Vec3({self.x:.12f}, {self.y:.12f}, {self.z:.12f})"

    def __getitem__(self, index: int) -> float:
        return (self.x, self.y, self.z)[index]

    def __setitem__(self, index: int, value: int | float) -> None:
        if not isinstance(value, (int, float)):
            raise TypeError("value must be an int or float")
        if index in (0, -3):
            self.x = value
        elif index in (1, -2):
            self.y = value
        elif index in (2, -1):
            self.z = value
        else:
            raise IndexError("index out of range")

    def __eq__(self, other: Vec3) -> bool:
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other: Vec3) -> bool:
        return not self == other

    def __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    def __add__(self, other: Vec3) -> Vec3:
        other_class_check_for_operand(self, other, "+")
        return self.__class__(self.x + other.x, self.y + other.y, self.z + other.z)

    def __radd__(self, other: Vec3) -> Vec3:
        return self.__add__(other)

    def __iadd__(self, other: Vec3) -> Vec3:
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __sub__(self, other: Vec3) -> Vec3:
        other_class_check_for_operand(self, other, "-")
        return self.__class__(self.x - other.x, self.y - other.y, self.z - other.z)

    def __rsub__(self, other: Vec3) -> Vec3:
        return self.__class__(other.x - self.x, other.y - self.y, other.z - self.z)

    def __isub__(self, other: Vec3) -> Vec3:
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __mul__(self, value: int | float | Vec3) -> Vec3:
        if not isinstance(value, (int, float, Vec3)):
            raise TypeError("value must be an int, float or Vec3")
        if isinstance(value, (int, float)):
            return Vec3(self.x * value, self.y * value, self.z * value)

        return Vec3(self.x * value.x, self.y * value.y, self.z * value.z)

    def __rmul__(self, value: int | float | Vec3) -> Vec3:
        return self.__mul__(value)

    def __imul__(self, value: int | float | Vec3) -> Vec3:
        if isinstance(value, (int, float)):
            self.x *= value
            self.y *= value
            self.z *= value
            return self
        if isinstance(value, Vec3):
            self.x *= value.x
            self.y *= value.y
            self.z *= value.z
            return self
        raise TypeError("value must be an int, float or Vec3")

    def __truediv__(self, value: int | float | Vec3) -> Vec3:
        if not isinstance(value, (int, float, Vec3)):
            raise TypeError("value must be an int, float or Vec3")
        if isinstance(value, (int, float)):
            return Vec3(self.x / value, self.y / value, self.z / value)

        return Vec3(self.x / value.x, self.y / value.y, self.z / value.z)

    def __itruediv__(self, value: int | float | Vec3) -> Vec3:
        if isinstance(value, (int, float)):
            self.x /= value
            self.y /= value
            self.z /= value
            return self
        if isinstance(value, Vec3):
            self.x /= value.x
            self.y /= value.y
            self.z /= value.z
            return self
        raise TypeError("value must be an int, float or Vec3")

    def __floordiv__(self, value: int | float | Vec3) -> Vec3:
        if not isinstance(value, (int, float, Vec3)):
            raise TypeError("value must be an int, float or Vec3")
        if isinstance(value, (int, float)):
            return Vec3(self.x // value, self.y // value, self.z // value)

        return Vec3(self.x // value.x, self.y // value.y, self.z // value.z)

    def __ifloordiv__(self, value: int | float | Vec3) -> Vec3:
        if isinstance(value, (int, float)):
            self.x //= value
            self.y //= value
            self.z //= value
            return self
        if isinstance(value, Vec3):
            self.x //= value.x
            self.y //= value.y
            self.z //= value.z
            return self
        raise TypeError("value must be an int, float or Vec3")

    def __mod__(self, value: int | float | Vec3) -> Vec3:
        if not isinstance(value, (int, float, Vec3)):
            raise TypeError("value must be an int, float or Vec3")
        if isinstance(value, (int, float)):
            return Vec3(self.x % value, self.y % value, self.z % value)

        return Vec3(self.x % value.x, self.y % value.y, self.z % value.z)

    def __imod__(self, value: int | float | Vec3) -> Vec3:
        if isinstance(value, (int, float)):
            self.x %= value
            self.y %= value
            self.z %= value
            return self
        if isinstance(value, Vec3):
            self.x %= value.x
            self.y %= value.y
            self.z %= value.z
            return self
        raise TypeError("value must be an int, float or Vec3")

    def __len__(self) -> int:
        return 3

    def copy(self) -> Vec3:
        return deepcopy(self)

    def dot(self, other: vec_type) -> float:
        other = vectorize_arg(other)
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: vec_type) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def distance_to(self, other: Vec3) -> float:
        other = vectorize_arg(other)
        return math.sqrt(
            (self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2
        )

    def isclose(self, other: Vec3, rel_tol: float = 1e-5, abs_tol: float = 0.0) -> bool:
        other = vectorize_arg(other)
        self_coords = [self.x, self.y, self.z]
        other_coords = [other.x, other.y, other.z]
        return all(
            math.isclose(
                self_coords[i],
                other_coords[i],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            )
            for i in range(3)
        )

    def norm(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normsq(self) -> float:
        return self.x**2 + self.y**2 + self.z**2

    def normalize(self) -> None:
        self /= self.norm()

    def angle(self, other: Vec3) -> float:
        if not isinstance(other, Vec3):
            raise TypeError(f"Expected Vec3, got {type(other)}")
        other = vectorize_arg(other)
        dot = self.dot(other)
        cos_theta = dot / (self.norm() * other.norm())
        return math.acos(cos_theta)

    def mid_point(self, other: Vec3) -> Vec3:
        return (self + other) / 2

    def translate(self, trans_vec) -> None:
        trans_vec = vectorize_arg(trans_vec)
        self.x += trans_vec.x
        self.y += trans_vec.y
        self.z += trans_vec.z

    def mirror(self, sx: int, sy: int, sz: int) -> None:
        if sx not in (-1, 1) or sy not in (-1, 1) or sz not in (-1, 1):
            raise ValueError("mirror factors must be 1 or -1")
        self.x *= sx
        self.y *= sy
        self.z *= sz

    def mirror_by_plane(self, p1: Vec3, p2: Vec3, p3: Vec3) -> None:
        if not all(isinstance(p, Vec3) for p in (p1, p2, p3)):
            raise TypeError("p1, p2, p3 must be Vec3 instances")

        v1 = p2 - p1
        v2 = p3 - p1

        normal = v1.cross(v2)

        if normal.norm() < 1e-6:
            raise ValueError("the three points are collinear")

        # calculate the projection of the point on the mirror plane
        projection = self - normal * (self - p1).dot(normal) / normal.normsq()

        # calculate the mirror image of the point
        self += 2 * (projection - self)

    def matmul(self, mat: mat_type) -> None:
        x_rot = mat[0][0] * self.x + mat[0][1] * self.y + mat[0][2] * self.z
        y_rot = mat[1][0] * self.x + mat[1][1] * self.y + mat[1][2] * self.z
        z_rot = mat[2][0] * self.x + mat[2][1] * self.y + mat[2][2] * self.z

        self.x = x_rot
        self.y = y_rot
        self.z = z_rot

    def rotate_by_axis(
        self,
        axis_point1: vec_type,
        axis_point2: vec_type,
        deg: float | int,
    ) -> None:
        """
        :param axis_point1: One point on the rotation axis
        :param axis_point2: Another point on the rotation axis
        :param deg: Rotation angle (degrees)
        """

        # if scipy and numpy is available, you can use the following code
        """
        angle_radians = np.deg2rad(deg)
        axis_vector = np.array(
            [
                axis_point2.x - axis_point1.x,
                axis_point2.y - axis_point1.y,
                axis_point2.z - axis_point1.z,
            ]
        )
        axis_unit_vector = axis_vector / np.linalg.norm(axis_vector)
        point = np.array([*self]) - np.array([*axis_point1])
        rotation = R.from_rotvec(angle_radians * axis_unit_vector)
        rotated_point = rotation.apply(point)
        rotated_point += np.array([*axis_point1])
        self.x, self.y, self.z = rotated_point.tolist()
        """

        axis_point1 = vectorize_arg(axis_point1)
        axis_point2 = vectorize_arg(axis_point2)

        angle_radians = math.radians(deg)
        axis_vector = axis_point2 - axis_point1
        axis_unit_vector = axis_vector / axis_vector.norm()

        # translate to move the rotation axis to the origin
        self -= axis_point1

        # calculate elements of rotation matrix
        cos_theta = math.cos(angle_radians)
        sin_theta = math.sin(angle_radians)
        one_minus_cos = 1 - cos_theta

        # Rodrigues' rotation formula
        x_rot = (
            (cos_theta + axis_unit_vector.x**2 * one_minus_cos) * self.x
            + (
                axis_unit_vector.x * axis_unit_vector.y * one_minus_cos
                - axis_unit_vector.z * sin_theta
            )
            * self.y
            + (
                axis_unit_vector.x * axis_unit_vector.z * one_minus_cos
                + axis_unit_vector.y * sin_theta
            )
            * self.z
        )

        y_rot = (
            (
                axis_unit_vector.y * axis_unit_vector.x * one_minus_cos
                + axis_unit_vector.z * sin_theta
            )
            * self.x
            + (cos_theta + axis_unit_vector.y**2 * one_minus_cos) * self.y
            + (
                axis_unit_vector.y * axis_unit_vector.z * one_minus_cos
                - axis_unit_vector.x * sin_theta
            )
            * self.z
        )

        z_rot = (
            (
                axis_unit_vector.z * axis_unit_vector.x * one_minus_cos
                - axis_unit_vector.y * sin_theta
            )
            * self.x
            + (
                axis_unit_vector.z * axis_unit_vector.y * one_minus_cos
                + axis_unit_vector.x * sin_theta
            )
            * self.y
            + (cos_theta + axis_unit_vector.z**2 * one_minus_cos) * self.z
        )

        # translate back to the original position
        self.x = x_rot + axis_point1.x
        self.y = y_rot + axis_point1.y
        self.z = z_rot + axis_point1.z


# can be used as a type hint. not isinstance() check
vec_type = Vec3 | list[int | float]
mat_type = list[Vec3 | list[int | float]]


def vectorize_arg(arg: vec_type) -> Vec3:
    if isinstance(arg, Vec3):
        return arg
    if isinstance(arg, list) and len(arg) == 3:
        if all(isinstance(i, (int, float)) for i in arg):
            return Vec3(*arg)
    raise ValueError("arg must be a Vec3 instance or a list of 3 int or float")


def other_class_check_for_operand(self, other, operand: str) -> None:
    if not isinstance(other, self.__class__):
        return NotImplementedError(
            f"unsupported operand type(s) for {operand}:"
            + f"'{self.__class__.__name__}' and '{other.__class__.__name__}'"
        )


def is_vec_type(vec: list[float]) -> bool:
    return (
        isinstance(vec, list)
        and len(vec) == 3
        and all(isinstance(i, (int, float)) for i in [*vec])
    )


def is_mat_type(mat: list[list[float]] | Mat3) -> bool:
    return isinstance(mat, Mat3) or (
        isinstance(mat, list)
        and len(mat) == 3
        and all(is_vec_type(vec) for vec in mat if isinstance(vec, list))
    )
