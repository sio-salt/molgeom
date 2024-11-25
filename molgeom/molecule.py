from __future__ import annotations
import yaml
import threading
import importlib.resources
import networkx as nx
import copy
from collections.abc import Iterable
from easyvec import Vec3
from molgeom.utils.fancy_indexing_list import FancyIndexingList
from molgeom.data.consts import ANGST2BOHR_GAU16, ATOMIC_NUMBER
from molgeom.atom import Atom
from molgeom.utils.decorators import args_to_set


_bond_data = None
_mole_data_lock = threading.Lock()


def _load_bond_data():
    global _bond_data
    if _bond_data is None:
        with _mole_data_lock:
            if _bond_data is None:
                with importlib.resources.open_text(
                    "molgeom.data", "bond_lengths.yaml"
                ) as f:
                    _bond_data = yaml.safe_load(f)
    return _bond_data


class Molecule:
    def __init__(self, *atoms: Atom | Iterable[Atom]) -> None:
        self.atoms: FancyIndexingList[Atom] = FancyIndexingList()
        self._data = _load_bond_data()
        if atoms:
            self.add_atoms(*atoms)

    def __len__(self) -> int:
        return len(self.atoms)

    def __eq__(self, other: Molecule) -> bool:
        return self.atoms == other.atoms

    def __ne__(self, other: Molecule) -> bool:
        return not self == other

    def __contains__(self, atom: Atom) -> bool:
        return atom in self.atoms

    def __getitem__(self, index: int | slice | Iterable) -> Atom | Molecule:
        result = self.atoms[index]
        if isinstance(result, Atom):
            return result
        elif isinstance(result, Iterable):
            return Molecule(*result)
        else:
            raise TypeError("Index must be int, slice, list, or tuple")

    def __setitem__(self, index: int, atom: Atom) -> None:
        self.atoms[index] = atom

    def __delitem__(self, idx: int) -> None:
        del self.atoms[idx]

    def __iter__(self):
        for atom in self.atoms:
            if not isinstance(atom, Atom):
                raise TypeError(
                    "Invalid element type: atom must be Atom object" + f"{type(atom)=}"
                )
            yield atom

    def __str__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    def __repr__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    def sort(self, key=None) -> None:
        if key is None:
            self.atoms.sort(key=lambda atom: ATOMIC_NUMBER[atom.symbol])
        else:
            self.atoms.sort(key=key)

    def pop(self, index: int) -> Atom:
        return self.atoms.pop(index)

    def copy(self) -> Molecule:
        return copy.deepcopy(self)

    @args_to_set
    def filter_by_symbols(self, symbols: str | Iterable[str]) -> Molecule:
        return Molecule(*[atom for atom in self if atom.symbol in symbols])

    def add_atoms(self, *atoms: Atom | Iterable[Atom]) -> None:
        """
        Add atoms to the molecule. Accepts either:
        - A list or tuple of Atom objects
        - Multiple Atom objects as separate arguments
        """
        if not atoms:
            return

        error_msg = "atoms must be a Atom object or a list/tuple of Atom objects"
        if isinstance(atoms, Atom):
            self.atoms.append(atoms)
            return
        elif isinstance(atoms, Iterable):
            for i, atom in enumerate(atoms):
                if not isinstance(atom, Atom):
                    raise TypeError(
                        f"invalid type: {i}th element type : {type(atom) = }\n{error_msg}"
                    )
            self.atoms.extend(atoms)
        else:
            raise TypeError(error_msg)

    def get_symbols(self) -> list[str]:
        return tuple(atom.symbol for atom in self)

    def get_formula(self) -> str:
        symbol_count = dict()
        for atom in self:
            symbol_count[atom.symbol] = symbol_count.get(atom.symbol, 0) + 1

        formula = "-".join(
            f"{symbol}{count}" if count > 1 else symbol
            for symbol, count in sorted(
                symbol_count.items(), key=lambda x: ATOMIC_NUMBER[x[0]]
            )
        )
        return formula

    def get_bonds(
        self,
        tol: float = 0.15,
    ) -> list[tuple[int, int]]:
        bonds = list()
        num_atoms = len(self)
        for i in range(num_atoms):
            ai = self[i]
            for j in range(i + 1, num_atoms):
                aj = self[j]
                if ai.is_bonded_to(aj, tol):
                    bonds.append((i, j))
        return bonds

    def get_bond_clusters(self, tol: int = 0.15) -> list[Molecule]:
        G = nx.Graph()
        G.add_nodes_from(range(len(self)))
        G.add_edges_from(self.get_bonds(tol))
        copied_mol = self.copy()
        return [copied_mol[list(cluster)] for cluster in nx.connected_components(G)]

    def get_cycles(self, length_bound: int = None, tol: float = 0.15) -> list[Molecule]:
        G = nx.Graph()
        G.add_edges_from(self.get_bonds(tol))
        cycles = nx.simple_cycles(G, length_bound=length_bound)
        return [self[list(cycle)] for cycle in cycles]

    def translate(self, trans_vec: Vec3) -> None:
        for atom in self.atoms:
            atom.translate(trans_vec)

    def mirror(self, sx: int, sy: int, sz: int) -> None:
        for atom in self.atoms:
            atom.mirror(sx, sy, sz)

    def mirror_by_plane(self, p1: Vec3, p2: Vec3, p3: Vec3) -> None:
        for atom in self.atoms:
            atom.mirror_by_plane(p1, p2, p3)

    def rotate_by_axis(
        self, axis_point1: Vec3, axis_point2: Vec3, angle_degrees: float
    ) -> None:
        """
        :param axis_point1: One point on the rotation axis
        :param axis_point2: Another point on the rotation axis
        :param angle_degrees: Rotation angle (degrees)
        """

        for atom in self.atoms:
            atom.rotate_by_axis(axis_point1, axis_point2, angle_degrees)

    def total_mass(self) -> float:
        return sum(atom.mass for atom in self)

    def center_of_mass(self) -> Vec3:
        total_mass = self.total_mass()
        com_x = sum(atom.x * atom.mass for atom in self) / total_mass
        com_y = sum(atom.y * atom.mass for atom in self) / total_mass
        com_z = sum(atom.z * atom.mass for atom in self) / total_mass
        return Vec3(com_x, com_y, com_z)

    def nuclear_repulsion(self) -> float:
        nuclrep = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                dist_angst = self.atoms[i].distance_to(self.atoms[j])
                dist_bohr = dist_angst * ANGST2BOHR_GAU16
                nuclrep += (
                    self.atoms[i].atomic_number
                    * self.atoms[j].atomic_number
                    / dist_bohr
                )
        return nuclrep

    def to_xyz(self) -> str:
        return "\n".join([atom.to_xyz() for atom in self])

    def write_to_xyz(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write(f"{len(self)}\n")
            f.write(self.get_formula() + "\n")
            f.write(self.to_xyz())
