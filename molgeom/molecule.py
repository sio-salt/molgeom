from __future__ import annotations
import copy
import networkx as nx
from collections.abc import Iterable
from easyvec import Vec3
from molgeom.utils.fancy_indexing_list import FancyIndexingList
from molgeom.data.consts import ANGST2BOHR_GAU16, ATOMIC_NUMBER
from molgeom.atom import Atom
from molgeom.utils.decorators import args_to_set, args_to_list


class Molecule:
    """
    A class to represent a molecule.
    """

    def __init__(
        self, *atoms: Atom, lattice_vecs: list[Vec3 | list[float | int]] | None = None
    ):
        self.atoms: FancyIndexingList[Atom] = FancyIndexingList()
        if atoms:
            self.add_atoms_from(atoms)
        self._lattice_vecs: list[Vec3] | None = None
        self.lattice_vecs = lattice_vecs
        self._bonds: list[tuple[int, int]] | None = None
        self._cycles: list[Molecule] | None = None

    @property
    def lattice_vecs(self) -> list[Vec3] | None:
        return self._lattice_vecs

    @lattice_vecs.setter
    def lattice_vecs(self, lattice_vecs: list[Vec3 | list[float | int]] | None) -> None:
        if lattice_vecs is None:
            self._lattice_vecs = None
        elif not all(isinstance(vec, (Vec3, list)) for vec in lattice_vecs):
            raise TypeError(
                "All elements must be Vec3 objects or list of 3 floats or ints"
            )
        elif not all(len(vec) == 3 for vec in lattice_vecs):
            raise ValueError("All lattice vectors must be of length 3")
        else:
            self._lattice_vecs = [
                Vec3(*vec) if isinstance(vec, list) else vec for vec in lattice_vecs
            ]

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

    @classmethod
    def sorted(cls, mols: Molecule | Iterable[Molecule], key=None) -> Molecule:
        """
        Create a new molecule by sorting the atoms of the input molecule(s).
        """
        if not all(isinstance(mol, Molecule) for mol in mols):
            raise TypeError("All elements must be Molecule objects")
        sorted_mols = cls()
        for mol in mols:
            mol.sort(key=key)
            sorted_mols.add_atoms_from(mol)

        return sorted_mols

    def sort(self, key=None) -> None:
        """
        Sort the atoms of the molecule.
        """
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

    @classmethod
    @args_to_list
    def from_atoms(cls, atoms: Iterable[Atom]) -> Molecule:
        """
        Create a new molecule from Atom objects.
        """
        if not all(isinstance(atom, Atom) for atom in atoms):
            raise TypeError("All elements must be Atom objects")
        return cls(*atoms)

    def add_atom(self, atom: Atom) -> None:
        self.atoms.append(atom)

    def add_atoms_from(self, atoms: Iterable[Atom]) -> None:
        """
        Add atoms to the molecule.
        Accepts Iterable containing one or more Atom objects.
        """
        if not isinstance(atoms, Iterable):
            raise TypeError("atoms must be an Iterable of Atom objects")
        not_atoms = [atom for atom in atoms if not isinstance(atom, Atom)]
        if not_atoms:
            raise TypeError(
                "Invalid element type: atoms must be Atom objects " + f"{not_atoms=}"
            )
        self.atoms.extend(atoms)

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

        if self._bonds is not None:
            return self._bonds

        bonds = list()
        num_atoms = len(self)

        for i in range(num_atoms):
            ai = self[i]
            for j in range(i + 1, num_atoms):
                aj = self[j]
                if ai.is_bonded_to(aj, tol):
                    bonds.append((i, j))
        self._bonds = bonds

        return bonds

    def get_clusters(self, tol: int = 0.15) -> list[Molecule]:
        G = nx.Graph()
        G.add_nodes_from(range(len(self)))
        G.add_edges_from(self.get_bonds(tol))
        copied_mol = self.copy()
        return [copied_mol[list(cluster)] for cluster in nx.connected_components(G)]

    def get_cycles(self, length_bound: int = None, tol: float = 0.15) -> list[Molecule]:
        if self._cycles is not None:
            return self._cycles
        G = nx.Graph()
        G.add_edges_from(self.get_bonds(tol))
        cycles = nx.simple_cycles(G, length_bound=length_bound)
        self._cycles = [self[list(cycle)] for cycle in cycles]
        return self._cycles

    @classmethod
    @args_to_list
    def merged(cls, mols: Molecule | Iterable[Molecule]) -> Molecule:
        """
        Create a new molecule by merging multiple Molecule objects.
        lattice vectors, bonds, and cycles are initialized to None
        """
        if not all(isinstance(mol, Molecule) for mol in mols):
            raise TypeError("All elements must be Molecule objects")
        merged = cls()
        for mol in mols:
            merged.add_atoms_from(mol)
        return merged

    @args_to_list
    def merge(self, mols: Molecule | Iterable[Molecule]) -> Molecule:
        """
        Merge other Molecule objects into this Molecule object.
        resets the lattice vectors, bonds, and cycles
        """
        if not all(isinstance(mol, Molecule) for mol in mols):
            raise TypeError(
                "All elements must be Molecule objects\n"
                + f"{[(i, type(mol)) for i, mol in enumerate(mols) if not isinstance(mol, Molecule)]}"
            )

        for mol in mols:
            self.add_atoms_from(mol)

        self._lattice_vecs = None
        self._bonds = None
        self._cycles = None

    def translate(self, trans_vec: Vec3 | list[float | int]) -> None:
        if not isinstance(trans_vec, (Vec3, list)):
            raise TypeError("trans_vec must be Vec3 object or list of 3 floats or ints")
        if isinstance(trans_vec, list):
            if len(trans_vec) != 3:
                raise ValueError("trans_vec must be of length 3")
            if not all(isinstance(i, (int, float)) for i in trans_vec):
                raise TypeError("elements of trans_vec must be int or float")
            trans_vec = Vec3(*trans_vec)

        for atom in self:
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

    def replicate(self, rep_a: list[int], rep_b: list[int], rep_c: list[int]) -> None:
        """
        Replicate the molecule in the a, b, and c directions.
        """
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to replicate the molecule.")

        if not all(isinstance(rep, list) for rep in (rep_a, rep_b, rep_c)):
            raise TypeError(
                "Replication vectors must be of length 2 list\n"
                + f"{(rep_a, rep_b, rep_c)=}"
            )

        if len(rep_a) != 2 or len(rep_b) != 2 or len(rep_c) != 2:
            raise ValueError(
                "Replication vectors must be of length 2 list\n"
                + f"{(rep_a, rep_b, rep_c)=}"
            )

        if rep_a[0] >= rep_a[1] or rep_b[0] >= rep_b[1] or rep_c[0] >= rep_c[1]:
            raise ValueError(
                "Replication vectors must be of the form [start, end], start < end\n"
                + f"{(rep_a, rep_b, rep_c)=}"
            )

        tmp_mol = self.copy()
        self.atoms = FancyIndexingList()
        for i in range(rep_a[0], rep_a[1]):
            for j in range(rep_b[0], rep_b[1]):
                for k in range(rep_c[0], rep_c[1]):

                    trans_vec = (
                        i * tmp_mol.lattice_vecs[0]
                        + j * tmp_mol.lattice_vecs[1]
                        + k * tmp_mol.lattice_vecs[2]
                    )
                    mol_copied = tmp_mol.copy()
                    mol_copied.translate(trans_vec)
                    self.merge(mol_copied)

    def total_mass(self) -> float:
        return sum(atom.mass for atom in self)

    def center_of_mass(self) -> Vec3:
        total_mass = self.total_mass()
        com_x = sum(atom.x * atom.mass for atom in self) / total_mass
        com_y = sum(atom.y * atom.mass for atom in self) / total_mass
        com_z = sum(atom.z * atom.mass for atom in self) / total_mass
        return Vec3(com_x, com_y, com_z)

    def nuclear_repulsion(self) -> float:
        """
        Calculate the nuclear repulsion energy of the molecule.
        uses gaussian16 constants
        """
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

    def nuclrep(self) -> float:
        """
        Alias for nuclear_repulsion method
        """
        return self.nuclear_repulsion()

    def electrostatic_energy(self) -> float:
        """
        Calculate the electrostatic energy of the molecule.
        """
        elec_energy = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                dist_angst = self.atoms[i].distance_to(self.atoms[j])
                dist_bohr = dist_angst * ANGST2BOHR_GAU16
                if self.atoms[i].charge is None or self.atoms[j].charge is None:
                    raise ValueError("All atoms must have their charge set.")
                elec_energy += self.atoms[i].charge * self.atoms[j].charge / dist_bohr
        return elec_energy

    def to_xyz(self) -> str:
        return "\n".join([atom.to_xyz() for atom in self])

    def to_dict(self) -> dict:
        mol_dict = dict()
        for i, atom in enumerate(self):
            mol_dict[i] = atom.to_dict()
        return mol_dict

    def write_to_xyz(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write(f"{len(self)}\n")
            f.write(self.get_formula() + "\n")
            f.write(self.to_xyz())
        print(f"File written to {filepath}")
