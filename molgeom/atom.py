from __future__ import annotations
from easyvec import Vec3
from .consts import ATOMIC_NUMBER, ATOMIC_MASSES


class Atom(Vec3):
    def __init__(self, symbol: str, x: float, y: float, z: float) -> None:
        super().__init__(x, y, z)
        self.symbol = symbol
        self.mass = ATOMIC_MASSES.get(self.symbol, 0.0)
        self.atomic_number = ATOMIC_NUMBER.get(self.symbol, 0)

    def __str__(self) -> str:
        return f"{self.symbol:2s} {self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    def __repr__(self) -> str:
        return f"Atom({self.symbol!r}, {self.x:.12f}, {self.y:.12f}, {self.z:.12f})"

    def __eq__(self, other: Atom) -> bool:
        return self.symbol == other.symbol and super().__eq__(other)

    def __ne__(self, other: Atom) -> bool:
        return not self == other

    def __lt__(self, other: Atom) -> bool:
        if self.atomic_number != other.atomic_number:
            return self.atomic_number < other.atomic_number
        if self.mass != other.mass:
            return self.mass < other.mass

    @classmethod
    def from_point(cls, symbol: str, point: Vec3) -> Atom:
        return cls(symbol, point.x, point.y, point.z)

    def to_xyz(self) -> str:
        return f"{self.symbol:2s} {self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    def closest_atom(self, mole, symbol=None) -> Atom:
        closest = mole.atoms[0]
        for atom in mole:
            if symbol is not None and atom.symbol != symbol:
                continue
            if self.distance_to(atom) < self.distance_to(closest):
                closest = atom
        return closest

