from __future__ import annotations
import yaml
import threading
import importlib.resources
from easyvec import Vec3


_pt_data = None
_atomic_data_lock = threading.Lock()


def _load_atomic_data():
    global _pt_data
    global _bond_data
    if _pt_data is None:
        with _atomic_data_lock:
            if _pt_data is None:
                with importlib.resources.open_text(
                    "molgeom.data", "periodic_table.yaml"
                ) as f:
                    _pt_data = yaml.load(f)
                with importlib.resources.open_text(
                    "molgeom.data", "bonds_jmol_ob.yaml"
                ) as f:
                    _bond_data = yaml.load(f)
    return _pt_data, _bond_data


class Atom(Vec3):
    def __init__(self, symbol: str, x: float, y: float, z: float) -> None:
        super().__init__(x, y, z)
        self.symbol = symbol
        self._data, self._std_bond_rad = self.get_atomic_data(symbol)
        self.mass = self._data.get("Atomic mass", 0.0)
        self.atomic_number = self._data.get("Atomic no", 0.0)

    @classmethod
    def from_point(cls, symbol: str, point: Vec3) -> Atom:
        return cls(symbol, point.x, point.y, point.z)

    @staticmethod
    def get_atomic_data(symbol):
        _pt_data, _bond_data = _load_atomic_data()
        return _pt_data.get(symbol, {}), _bond_data.get(symbol, {})

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

    def to_xyz(self) -> str:
        return f"{self.symbol:2s} {self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    def is_bonded_to(self, other, tol=0.15):
        dist_angst = self.distance_to(other)
        std_bond_len = self._std_bond_rad + other._std_bond_rad
        return abs(dist_angst - std_bond_len) <= tol

    def closest_atom(self, mole, symbol=None) -> Atom:
        closest = mole.atoms[0]
        for atom in mole:
            if symbol is not None and atom.symbol != symbol:
                continue
            if self.distance_to(atom) < self.distance_to(closest):
                closest = atom
        return closest
