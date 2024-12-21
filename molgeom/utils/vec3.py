from __future__ import annotations

from typing import Self, TypeVar

import numpy as np
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation


# can be used as a type hint. not isinstance() check
# Tvec is a type variable that can be bound to Vec3 or a subclass of Vec3
Tvec = TypeVar("Tvec", bound="Vec3")
scalar_type = int | float | np.integer | np.floating
scalar_type_tuple = (int, float, np.integer, np.floating)
VecLike = Tvec | list[scalar_type, scalar_type, scalar_type] | ArrayLike


class Vec3:
    def __init__(self, x: scalar_type, y: scalar_type, z: scalar_type) -> None:
        if not all(isinstance(i, scalar_type_tuple) for i in (x, y, z)):
            raise TypeError("x, y, z must be an int or float")

        self._coord = np.asarray([x, y, z], dtype=float)

    @property
    def coord(self) -> np.ndarray:
        return self._coord

    @coord.setter
    def coord(self, array: np.ndarray) -> None:
        if not isinstance(array, np.ndarray):
            raise TypeError("value must be a numpy.ndarray")
        if array.shape != (3,):
            raise ValueError("value must have shape (3,)")
        self._coord[:] = array

    @property
    def x(self) -> float:
        return self._coord[0]

    @x.setter
    def x(self, value: scalar_type) -> None:
        self._validate_scalar(value)
        self._coord[0] = float(value)

    @property
    def y(self) -> float:
        return self._coord[1]

    @y.setter
    def y(self, value: scalar_type) -> None:
        self._validate_scalar(value)
        self._coord[1] = float(value)

    @property
    def z(self) -> float:
        return self._coord[2]

    @z.setter
    def z(self, value: scalar_type) -> None:
        self._validate_scalar(value)
        self._coord[2] = float(value)

    def _validate_scalar(self, value: scalar_type) -> None:
        if not isinstance(value, scalar_type_tuple):
            raise TypeError(f"value must be an int or float, got {type(value)}")

    @classmethod
    def from_array(cls, array: ArrayLike) -> Vec3:
        array = np.asarray(array, dtype=float)
        if array.shape != (3,):
            raise ValueError(f"array must have shape (3,), got {array.shape}")
        return cls(*array)

    def to_list(self) -> list[float, float, float]:
        return [self.x, self.y, self.z]

    def to_dict(self) -> dict[str, float]:
        return {"x": self.x, "y": self.y, "z": self.z}

    def __str__(self) -> str:
        return f"{self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    def __repr__(self) -> str:
        return f"Vec3({self.x:.12f}, {self.y:.12f}, {self.z:.12f})"

    def __getitem__(self, index: int | slice) -> scalar_type | np.ndarray:
        return self.coord[index]

    def __setitem__(self, index: int | slice, value: scalar_type) -> None:
        self.coord[index] = value

    def __array__(self, dtype=None) -> np.ndarray:
        if dtype is not None:
            return self.coord.astype(dtype)
        return self._coord

    def __eq__(self, other: Self) -> bool:
        if not isinstance(other, type(self)):
            return False
        return np.array_equal(self.coord, other.coord)

    def __ne__(self, other: Self) -> bool:
        return not self == other

    def __neg__(self: Tvec) -> Tvec:
        return type(self).from_array(-self.coord)

    def __add__(self: Self, other: Tvec) -> Self:
        other_class_check_for_operand(self, other, "+")
        return type(self).from_array(self.coord + other.coord)

    def __radd__(self: Self, other: Tvec) -> Tvec:
        other_class_check_for_operand(self, other, "+")
        return type(other)(other.coord + self.coord)

    def __iadd__(self: Self, other: Tvec) -> Self:
        other_class_check_for_operand(self, other, "+")
        self.coord += other.coord
        return self

    def __sub__(self: Self, other: Tvec) -> Self:
        other_class_check_for_operand(self, other, "-")
        return type(self).from_array(self.coord - other.coord)

    def __rsub__(self: Self, other: Tvec) -> Tvec:
        other_class_check_for_operand(self, other, "-")
        return type(other)(other.coord - self.coord)

    def __isub__(self: Self, other: Tvec) -> Self:
        other_class_check_for_operand(self, other, "-")
        self.coord -= other.coord
        return self

    def __mul__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            return type(self).from_array(self.coord * value)
        if isinstance(value, self.__class__):
            return type(self).from_array(self.coord * value.coord)
        raise TypeError(f"value must be an scalar_type_tuple or Vec3{type(value)}")

    def __rmul__(self, value: scalar_type | Tvec) -> Tvec:
        return self.__mul__(value)

    def __imul__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            self.coord *= value
            return self
        if isinstance(value, self.__class__):
            self.coord *= value.coord
            return self
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __truediv__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            return type(self).from_array(self.coord / value)
        if isinstance(value, Vec3):
            return type(self).from_array(self.coord / value.coord)
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __rtruediv__(self, value: scalar_type | Tvec) -> Tvec:
        tmp = self.copy()
        tmp.coord[:] = value / tmp.coord
        return tmp

    def __itruediv__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            self.coord /= value
            return self
        if isinstance(value, Vec3):
            self.coord /= value.coord
            return self
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __floordiv__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            return type(self).from_array(self.coord // value)
        if isinstance(value, Vec3):
            return type(self).from_array(self.coord // value.coord)
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __rfloordiv__(self, value: scalar_type | Tvec) -> Tvec:
        tmp = self.copy()
        tmp.coord[:] = value // tmp.coord
        return tmp

    def __ifloordiv__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            self.coord //= value
            return self
        if isinstance(value, Vec3):
            self.coord //= value.coord
            return self
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __mod__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            return type(self).from_array(self.coord % value)
        if isinstance(value, Vec3):
            return type(self).from_array(self.coord % value.coord)
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __rmod__(self, value: scalar_type | Tvec) -> Tvec:
        tmp = self.copy()
        tmp.coord[:] = value % tmp.coord
        return tmp

    def __imod__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            self.coord %= value
            return self
        if isinstance(value, Vec3):
            self.coord %= value.coord
            return self
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __pow__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            return type(self).from_array(self.coord**value)
        if isinstance(value, Vec3):
            return type(self).from_array(self.coord**value.coord)
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __rpow__(self, value: scalar_type | Tvec) -> Tvec:
        tmp = self.copy()
        tmp.coord[:] = value**tmp.coord
        return tmp

    def __ipow__(self, value: scalar_type | Tvec) -> Self:
        if isinstance(value, scalar_type_tuple):
            self.coord **= value
            return self
        if isinstance(value, Vec3):
            self.coord **= value.coord
            return self
        raise TypeError(
            f"value must be an scalar_type_tuple or Vec3, got {type(value)}"
        )

    def __matmul__(self: Self, other: Tvec | np.ndarray) -> Tvec:
        if isinstance(other, Vec3):
            return type(self).from_array(self.coord @ other.coord)
        if isinstance(other, np.ndarray):
            return type(self).from_array(self.coord @ other)
        raise TypeError(f"value must be an Vec3 or np.ndarray, got {type(other)}")

    def __rmatmul__(self: Self, other: Tvec | np.ndarray) -> Tvec:
        if isinstance(other, Vec3):
            return type(other).from_array(other.coord @ self.coord)
        if isinstance(other, np.ndarray):
            return type(other).from_array(other @ self.coord)
        raise TypeError(f"value must be an Vec3 or np.ndarray, got {type(other)}")

    def __imatmul__(self: Self, other: Tvec | np.ndarray) -> Self:
        if isinstance(other, Vec3):
            self.coord @= other.coord
            return self
        if isinstance(other, np.ndarray):
            self.coord @= other
            return self
        raise TypeError(f"value must be an Vec3 or np.ndarray, got {type(other)}")

    def matmul(self, mat: np.ndarray) -> None:
        self.coord[:] = mat @ self.coord

    def __len__(self) -> int:
        return 3

    def copy(self) -> Self:
        return type(self).from_array(self._coord.copy())

    def dot(self, other: Tvec) -> float:
        return np.dot(self.coord, other.coord)

    def cross(self, other: Tvec) -> Self:
        return type(self).from_array(np.cross(self.coord, other.coord))

    def distance_to(self, other: VecLike) -> float:
        other = vec3fy(other)
        return np.linalg.norm(self.coord - other.coord)

    def distance_to_sq(self, other: VecLike) -> float:
        other = vec3fy(other)
        return np.sum((self.coord - other.coord) ** 2)

    def isclose(
        self, other: VecLike, rel_tol: float = 1e-5, abs_tol: float = 0.0
    ) -> bool:
        other = vec3fy(other)
        return np.allclose(self.coord, other.coord, rtol=rel_tol, atol=abs_tol)

    def norm(self) -> float:
        return np.linalg.norm(self.coord)

    def normsq(self) -> float:
        return np.sum(self.coord**2)

    def normalize(self) -> None:
        self /= self.norm()

    def normalized(self) -> Self:
        return self / self.norm()

    def angle_to(self, other: VecLike, deg: bool = True) -> float:
        other = vec3fy(other)
        angle = np.arccos(self.dot(other) / (self.norm() * other.norm()))
        return np.degrees(angle) if deg else angle

    def translate(self, trans: VecLike) -> None:
        self += vec3fy(trans)

    def scale(self, scale: scalar_type | VecLike) -> None:
        self *= scale

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

    def rotate_by_axis(
        self, axis_point1: VecLike, axis_point2: VecLike, deg: scalar_type
    ) -> None:
        axis_point1 = np.asarray(axis_point1)
        axis_point2 = np.asarray(axis_point2)
        angle_radians = np.deg2rad(deg)

        axis_vector = axis_point2 - axis_point1
        axis_unit_vector = axis_vector / axis_vector.norm()
        rot = Rotation.from_rotvec(angle_radians * axis_unit_vector)
        point = self.coord - axis_point1
        rotated_point = rot.apply(point)
        self._coord[:] = rotated_point + axis_point1


def vec3fy(arg: VecLike) -> Vec3:
    if isinstance(arg, Vec3):
        return arg
    if isinstance(arg, list) and len(arg) == 3:
        if all(isinstance(i, scalar_type_tuple) for i in arg):
            return Vec3(*arg)
    raise ValueError(
        f"arg must be a Vec3 instance or a list of 3 int or float, got {arg}"
    )


def other_class_check_for_operand(self, other, operand: str) -> None:
    # if other is not the same class as self or a subclass of self, raise an error
    if not isinstance(other, self.__class__):
        return NotImplementedError(
            f"unsupported operand type(s) for {operand}: \n"
            + f"{type(self).__name__} and {type(other).__name__}"
        )
