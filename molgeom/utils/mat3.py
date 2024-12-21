from __future__ import annotations
from copy import deepcopy

from .vec3 import Vec3


class Mat3:
    def __init__(self, mat: mat_type) -> None:

        if not is_mat_type(mat):
            raise TypeError(f"Expected 3x3 matrix list, got {type(mat)}")
        self.mat = [Vec3(*row) for row in mat]

    def __str__(self) -> str:
        return f"{str(self.mat[0])}\n{str(self.mat[1])}\n{str(self.mat[2])}"

    def __repr__(self) -> str:
        return f"Mat3({repr(self.mat[0])}\n{repr(self.mat[1])}\n{repr(self.mat[2])})"

    def __eq__(self, other: Mat3) -> bool:
        return self.mat == other.mat

    def __ne__(self, other: Mat3) -> bool:
        return not self == other

    def __getitem__(self, index: int) -> Vec3:
        return self.mat[index]

    def __setitem__(self, index: int, value: vec_type) -> None:
        if not is_vec_type(value):
            raise TypeError(f"Expected 3 len vector, got {type(value)}")
        if isinstance(value, Vec3):
            self.mat[index] = value
        else:
            self.mat[index] = Vec3(*value)

    def __iter__(self) -> Vec3:
        for i in range(3):
            yield Vec3(*self.mat[i])

    def __len__(self) -> int:
        return 3

    def __neg__(self) -> Mat3:
        return Mat3([[-self.mat[i][j] for j in range(3)] for i in range(3)])

    def __add__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        return Mat3(
            [[self.mat[i][j] + other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __radd__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        return Mat3(
            [[self.mat[i][j] + other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __iadd__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] += other.mat[i][j]
        return self

    def __sub__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        return Mat3(
            [[self.mat[i][j] - other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __rsub__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        return Mat3(
            [[self.mat[i][j] - other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __isub__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] -= other.mat[i][j]
        return self

    def __mul__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] * other for j in range(3)] for i in range(3)])

    def __rmul__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] * other for j in range(3)] for i in range(3)])

    def __imul__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] *= other
        return self

    def __matmul__(self, other: Mat3 | Vec3) -> Mat3 | Vec3:
        if not isinstance(other, (Mat3, Vec3)):
            raise TypeError(f"Expected Mat3 or Vec3, got {type(other)}")

        if isinstance(other, Vec3):
            return Vec3.from_array(
                [sum(self.mat[i][j] * other[j] for j in range(3)) for i in range(3)]
            )
        else:
            return Mat3(
                [
                    [
                        sum(self.mat[i][k] * other.mat[k][j] for k in range(3))
                        for j in range(3)
                    ]
                    for i in range(3)
                ]
            )

    def __rmatmul__(self, other: Mat3 | Vec3) -> Mat3 | Vec3:
        if not isinstance(other, (Mat3, Vec3)):
            raise TypeError(f"Expected Mat3 or Vec3, got {type(other)}")

        if isinstance(other, Vec3):
            return Vec3.from_array(
                [sum(self.mat[i][j] * other[j] for j in range(3)) for i in range(3)]
            )
        else:
            return Mat3(
                [
                    [
                        sum(self.mat[i][k] * other.mat[k][j] for k in range(3))
                        for j in range(3)
                    ]
                    for i in range(3)
                ]
            )

    def __imatmul__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise TypeError(f"Expected 3x3 matrix, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] = sum(self.mat[i][k] * other.mat[k][j] for k in range(3))
        return self

    def __truediv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] / other for j in range(3)] for i in range(3)])

    def __itruediv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] /= other
        return self

    def __floordiv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] // other for j in range(3)] for i in range(3)])

    def __ifloordiv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] //= other
        return self

    def __mod__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] % other for j in range(3)] for i in range(3)])

    def __imod__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] %= other
        return self

    def __pow__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        return Mat3([[self.mat[i][j] ** other for j in range(3)] for i in range(3)])

    def __ipow__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        for i in range(3):
            for j in range(3):
                self.mat[i][j] **= other
        return self

    def copy(self) -> Mat3:
        return Mat3(deepcopy(self.mat))

    def to_list(self) -> list[list[float]]:
        return deepcopy(self.mat)

    def matrix_power(self, other: int | float) -> Mat3:
        """calculate the power of the matrix, self**3 = self*self*self"""
        if not isinstance(other, (int, float)):
            raise TypeError(f"Expected scalar, got {type(other)}")
        if other == 0:
            return Mat3([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        if other < 0:
            raise ValueError("minus power not supported yet")
        tmp = self.copy()
        for _ in range(other - 1):
            tmp = tmp @ self
        return tmp

    def T(self) -> Mat3:
        return Mat3([[self.mat[j][i] for j in range(3)] for i in range(3)])

    def det(self) -> float:
        return (
            self[0][0] * self[1][1] * self[2][2]
            + self[0][1] * self[1][2] * self[2][0]
            + self[0][2] * self[1][0] * self[2][1]
            - self[0][2] * self[1][1] * self[2][0]
            - self[0][0] * self[1][2] * self[2][1]
            - self[0][1] * self[1][0] * self[2][2]
        )

    def inv(self) -> Mat3:
        """
        mat^(-1) = 1/det(mat) * [[a22a33 - a23a32, -(a12a33 - a13a32), a12a23 - a13a22],
                                 [-(a21a33 - a23a31), a11a33 - a13a31, -(a11a23 - a13a21)],
                                 [a21a32 - a22a31, -(a11a32 - a12a31), a11a22 - a12a21]]
        """
        det = self.det()
        if det == 0:
            raise ValueError("Matrix is not invertible")
        return Mat3(
            [
                [
                    (self[1][1] * self[2][2] - self[1][2] * self[2][1]) / det,
                    -(self[0][1] * self[2][2] - self[0][2] * self[2][1]) / det,
                    (self[0][1] * self[1][2] - self[0][2] * self[1][1]) / det,
                ],
                [
                    -(self[1][0] * self[2][2] - self[1][2] * self[2][0]) / det,
                    (self[0][0] * self[2][2] - self[0][2] * self[2][0]) / det,
                    -(self[0][0] * self[1][2] - self[0][2] * self[1][0]) / det,
                ],
                [
                    (self[1][0] * self[2][1] - self[1][1] * self[2][0]) / det,
                    -(self[0][0] * self[2][1] - self[0][1] * self[2][0]) / det,
                    (self[0][0] * self[1][1] - self[0][1] * self[1][0]) / det,
                ],
            ]
        )


vec_type = Vec3 | list[int | float]
mat_type = list[Vec3 | list[int | float]]


def is_vec_type(vec: list[float]) -> bool:
    return (
        isinstance(vec, (list, Vec3))
        and len(vec) == 3
        and all(isinstance(i, (int, float)) for i in [*vec])
    )


def is_mat_type(mat: list[list[float]] | Mat3) -> bool:
    return isinstance(mat, Mat3) or (
        isinstance(mat, list)
        and len(mat) == 3
        and all(is_vec_type(vec) for vec in mat if isinstance(vec, list))
    )
