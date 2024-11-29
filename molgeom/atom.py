from __future__ import annotations
import json
import yaml
import threading
import importlib.resources
from typing import Any
from easyvec import Vec3
from molgeom.data import consts


_pt_data = None
_bond_pair_data = None
_atomic_data_lock = threading.Lock()


def _load_atomic_data():
    """
    Load atomic data from yaml files to _pt_data
    (periodic table data) and _bond_rad_data (bond data)
    This function is thread safe and will only load the data once
    """
    global _pt_data
    global _bond_rad_data
    if _pt_data is None:
        with _atomic_data_lock:
            if _pt_data is None:
                with importlib.resources.open_text(
                    "molgeom.data", "periodic_table.json"
                ) as f:
                    _pt_data = json.load(f)
                with importlib.resources.open_text(
                    "molgeom.data", "bonds_jmol_ob.yaml"
                ) as f:
                    _bond_rad_data = yaml.safe_load(f)
    return _pt_data, _bond_rad_data


def _load_bond_pair_data():
    """
    Load bond pair data from json file to _bond_pair_data
    This function is thread safe and will only load the data once
    """
    global _bond_pair_data
    if _bond_pair_data is None:
        with _atomic_data_lock:
            if _bond_pair_data is None:
                with importlib.resources.open_text(
                    "molgeom.data", "bond_pairs_len.json"
                ) as f:
                    _bond_pair_data = json.load(f)
    return _bond_pair_data


class Atom(Vec3):
    def __init__(self, symbol: str, x: float, y: float, z: float) -> None:
        super().__init__(x, y, z)
        if symbol not in consts.ATOMIC_NUMBER:
            raise ValueError(f"Invalid atomic symbol: {symbol}")
        self.symbol = symbol
        self._data, self._std_bond_rad = self.get_atomic_data(self.symbol)
        self._bond_pairs = _load_bond_pair_data()
        self.mass = self._data.get("Atomic mass", 0.0)
        self.atomic_number = self._data.get("Atomic no", 0.0)

    @classmethod
    def from_vec(cls, symbol: str, point: list | tuple | Vec3) -> Atom:
        if (
            isinstance(point, (list, tuple))
            and len(point) == 3
            and all(isinstance(i, (int, float)) for i in point)
        ):
            return cls(symbol, point[0], point[1], point[2])
        if isinstance(point, Vec3):
            return cls(symbol, point.x, point.y, point.z)
        return cls(symbol, point.x, point.y, point.z)

    @staticmethod
    def get_atomic_data(symbol):
        _pt_data, _bond_rad_data = _load_atomic_data()
        return _pt_data.get(symbol, {}), _bond_rad_data.get(symbol, {})

    def __getitem__(self, index: int) -> Any:
        return (self.symbol, self.x, self.y, self.z)[index]

    def __eq__(self, other: Atom) -> bool:
        return self.symbol == other.symbol and super().__eq__(other)

    def __ne__(self, other: Atom) -> bool:
        return not self == other

    def __lt__(self, other: Atom) -> bool:
        if self.atomic_number != other.atomic_number:
            return self.atomic_number < other.atomic_number
        if self.mass != other.mass:
            return self.mass < other.mass

    def __str__(self) -> str:
        return (
            f"Atom({self.symbol:2s}, {self.x:19.12f}, {self.y:19.12f}, {self.z:19.12f})"
        )

    def __repr__(self) -> str:
        return f"Atom({self.symbol!r}, {self.x:.12f}, {self.y:.12f}, {self.z:.12f})"

    def to_xyz(self) -> str:
        return f"{self.symbol:2s} {self.x:19.12f} {self.y:19.12f} {self.z:19.12f}"

    # def is_bonded_to(self, other, tol=0.12, lower_bound=None, upper_bound=None) -> bool:
    def is_bonded_to(self, other, tol=0.15) -> bool:
        dist_angst = self.distance_to(other)
        _bond_pair_data = _load_bond_pair_data()
        estimated_bond_len = self._std_bond_rad + other._std_bond_rad
        std_bond_lens = _bond_pair_data.get(self.symbol, {}).get(
            other.symbol, [estimated_bond_len]
        )

        # if lower_bound is not None and upper_bound is not None:
        #     for std_bond_len in std_bond_lens:
        #         if lower_bound <= std_bond_len <= upper_bound:
        #             return True
        for std_bond_len in std_bond_lens:
            if dist_angst <= std_bond_len + tol:
                return True
        return False

    def closest_atom(self, mole, symbol=None) -> Atom:
        closest = mole.atoms[0]
        for atom in mole:
            if symbol is not None and atom.symbol != symbol:
                continue
            if self.distance_to(atom) < self.distance_to(closest):
                closest = atom
        return closest
