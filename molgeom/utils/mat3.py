from __future__ import annotations
from copy import deepcopy

from .vec3 import Vec3


class Mat3:
    def __init__(self, mat: list[list[float]]):
        if not is_mat_type(mat):
            raise ValueError("Matrix must be 3x3")
        self.mat = mat

    def __str__(self) -> str:
        return f"{self.mat}"

    def __repr__(self) -> str:
        return f"Mat3({self.mat})"

    def __eq__(self, other: Mat3) -> bool:
        return self.mat == other.mat

    def __ne__(self, other: Mat3) -> bool:
        return not self == other

    def __getitem__(self, key: int) -> Vec3:
        return Vec3(*self.mat[key])

    def __setitem__(self, key: int, value: list[float]) -> None:
        if not is_vec_type(value):
            raise ValueError("Invalid vector type")
        self.mat[key] = value

    def __iter__(self) -> Vec3:
        for i in range(3):
            yield Vec3(*self.mat[i])

    def __add__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise ValueError("Invalid matrix type")
        return Mat3(
            [[self.mat[i][j] + other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __iadd__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise ValueError("Invalid matrix type")
        self.mat = [
            [self.mat[i][j] + other.mat[i][j] for j in range(3)] for i in range(3)
        ]
        return self

    def __sub__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise ValueError("Invalid matrix type")
        return Mat3(
            [[self.mat[i][j] - other.mat[i][j] for j in range(3)] for i in range(3)]
        )

    def __isub__(self, other: Mat3) -> Mat3:
        if not isinstance(other, Mat3):
            raise ValueError("Invalid matrix type")
        self.mat = [
            [self.mat[i][j] - other.mat[i][j] for j in range(3)] for i in range(3)
        ]
        return self

    def __mul__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        return Mat3([[self.mat[i][j] * other for j in range(3)] for i in range(3)])

    def __imul__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        self.mat = [[self.mat[i][j] * other for j in range(3)] for i in range(3)]
        return self

    def __matmul__(self, other: Mat3 | Vec3) -> Mat3 | Vec3:
        if not isinstance(other, (Mat3, Vec3)):
            raise TypeError("argument must be Mat3 or Vec3")

        if isinstance(other, Vec3):
            return Vec3.from_list(
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

    def __truediv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        return Mat3([[self.mat[i][j] / other for j in range(3)] for i in range(3)])

    def __itruediv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        self.mat = [[self.mat[i][j] / other for j in range(3)] for i in range(3)]
        return self

    def __floordiv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        return Mat3([[self.mat[i][j] // other for j in range(3)] for i in range(3)])

    def __ifloordiv__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        self.mat = [[self.mat[i][j] // other for j in range(3)] for i in range(3)]
        return self

    def __mod__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        return Mat3([[self.mat[i][j] % other for j in range(3)] for i in range(3)])

    def __imod__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        self.mat = [[self.mat[i][j] % other for j in range(3)] for i in range(3)]
        return self

    def __pow__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        return Mat3([[self.mat[i][j] ** other for j in range(3)] for i in range(3)])

    def __ipow__(self, other: int | float) -> Mat3:
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
        self.mat = [[self.mat[i][j] ** other for j in range(3)] for i in range(3)]
        return self

    def copy(self) -> Mat3:
        return Mat3(deepcopy(self.mat))

    def to_list(self) -> list[list[float]]:
        return deepcopy(self.mat)

    def matrix_power(self, other: int | float) -> Mat3:
        """calculate the power of the matrix, self**3 = self*self*self"""
        if not isinstance(other, (int, float)):
            raise ValueError("Invalid scalar type")
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
                    (self[1][1] * self[2][2] - self[1][2] * self[2][1]),
                    -(self[0][1] * self[2][2] - self[0][2] * self[2][1]),
                    (self[0][1] * self[1][2] - self[0][2] * self[1][1]),
                ],
                [
                    -(self[1][0] * self[2][2] - self[1][2] * self[2][0]),
                    (self[0][0] * self[2][2] - self[0][2] * self[2][0]),
                    -(self[0][0] * self[1][2] - self[0][2] * self[1][0]),
                ],
                [
                    (self[1][0] * self[2][1] - self[1][1] * self[2][0]),
                    -(self[0][0] * self[2][1] - self[0][1] * self[2][0]),
                    (self[0][0] * self[1][1] - self[0][1] * self[1][0]),
                ],
            ]
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
