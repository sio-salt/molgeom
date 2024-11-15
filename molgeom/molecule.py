from __future__ import annotations
from easyvec import Vec3
from molgeom._fancy_indexing_list import _FancyIndexingList
from molgeom.data.consts import ANGST2BOHR_GAU16
from molgeom.atom import Atom


class Molecule:
    def __init__(self, *atoms) -> None:
        self.atoms: _FancyIndexingList[Atom] = _FancyIndexingList()
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

    def __getitem__(self, index: int) -> Atom:
        return self.atoms[index]

    def __setitem__(self, index: int, atom: Atom) -> None:
        self.atoms[index] = atom

    def __iter__(self):
        return iter(self.atoms)

    def __str__(self):
        return f"class Molecule containing {len(self)} Atoms"

    def __repr__(self):
        return f"class Molecule containing {len(self)} Atoms"

    def sort(self) -> None:
        self.atoms.sort()

    def get_atoms_by_symbol(self, symbol: str) -> Molecule:
        matching_atoms = [atom for atom in self if atom.symbol == symbol]
        return Molecule(*matching_atoms)

    def add_atoms(self, *atoms) -> None:
        if not atoms:
            raise ValueError("atoms must not be empty")
        if not all(isinstance(atom, Atom) for atom in atoms):
            raise TypeError("atoms must be a list of Atom objects")
        for atom in atoms:
            self.atoms.append(atom)

    def append(self, *atoms) -> None:
        self.add_atoms(*atoms)

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

    def get_bonds(self, tol=0.15) -> list[tuple[int, int]]:
        bonds = list()
        num_atoms = len(self)
        for i in range(num_atoms):
            ai = self[i]
            for j in range(i + 1, num_atoms):
                aj = self[j]
                if ai.is_bonded_to(aj):
                    bonds.append((i, j))
        return bonds

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
