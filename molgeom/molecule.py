from __future__ import annotations

import copy
from collections.abc import Iterable

import networkx as nx
import numpy as np
from numpy.typing import ArrayLike
from cachetools import cachedmethod, LRUCache
from cachetools.keys import hashkey
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from molgeom.data.consts import ATOM_SYMBOLS
from molgeom.utils.fancy_indexing_list import FancyIndexingList
from molgeom.utils.vec3 import Vec3, VecLike
from molgeom.utils.decorators import args_to_list, args_to_set
from molgeom.utils.lattice_utils import cart2frac, frac2cart, lat_vecs_to_lat_params
from molgeom.utils.symmetry_utils import symmop_from_xyz_str
from molgeom.atom import Atom


default_tol = 0.15


class Molecule:
    """
    A class to represent a molecule.
    """

    def __init__(self, *atoms: Atom, lattice_vecs: ArrayLike | None = None):
        self.atoms: FancyIndexingList[Atom] = FancyIndexingList()
        self._symbols: np.ndarray = np.zeros((0,), dtype="<U2")
        self._coords: np.ndarray = np.zeros((0, 3), dtype=np.float64)
        if atoms:
            self.add_atoms_from(atoms)

        self._lattice_vecs: np.ndarray | None = None
        self.lattice_vecs: np.ndarray | None = lattice_vecs
        self._cache = LRUCache(maxsize=10)

    @property
    def symbols(self) -> np.ndarray:
        return self._symbols

    @symbols.setter
    def symbols(self, symbols: np.ndarray) -> None:
        if not isinstance(symbols, np.ndarray):
            raise TypeError("symbols must be a numpy array")
        self._symbols = symbols

    @property
    def coords(self) -> np.ndarray:
        return self._coords

    @coords.setter
    def coords(self, coords: np.ndarray) -> None:
        if not isinstance(coords, np.ndarray):
            raise TypeError("coords must be a numpy array")
        if coords.shape[1] != 3:
            raise ValueError(f"coords must be a Nx3 matrix, got shape {coords.shape}")
        self._coords = coords

    @property
    def lattice_vecs(self) -> np.ndarray | None:
        return self._lattice_vecs

    @lattice_vecs.setter
    def lattice_vecs(self, lattice_vecs: ArrayLike | None) -> None:
        if lattice_vecs is None:
            self._lattice_vecs = None
            return

        lat_vecs = np.asarray(lattice_vecs)
        if lat_vecs.shape != (3, 3):
            raise ValueError(
                f"lattice vectors must be a 3x3 matrix, got shape {lat_vecs.shape}"
            )
        if np.linalg.matrix_rank(lat_vecs) != 3:
            raise ValueError("lattice vectors must be linearly independent (rank 3)")
        if np.linalg.det(lat_vecs) <= 0:
            raise ValueError(
                f"lattice_vecs must form a right-handed coordinate system (det > 0), got {np.linalg.det(lat_vecs)}"
            )

        self._lattice_vecs = lat_vecs

    def __len__(self) -> int:
        return len(self.atoms)

    def __eq__(self, other: Molecule) -> bool:
        return self.atoms == other.atoms

    def __ne__(self, other: Molecule) -> bool:
        return not self == other

    def __contains__(self, atom: Atom) -> bool:
        if not isinstance(atom, Atom):
            raise TypeError("atom must be an Atom object")
        return atom in self.atoms

    def __getitem__(self, index: int | slice | Iterable) -> Atom | Molecule:
        result = self.atoms[index]
        if isinstance(result, Atom):
            return result
        elif isinstance(result, Iterable):
            return Molecule(*result)
        else:
            raise TypeError("Index must be int, slice, list, or tuple")

    # def __setitem__(self, key: int | slice, atoms: Atom | Molecule) -> None:
    #     if not isinstance(atoms, (Atom, Molecule)):
    #         raise TypeError(
    #             f"atoms must be an Atom or Molecule object, got {type(atoms).__name__}"
    #         )
    #     if isinstance(atoms, Molecule):
    #         atoms = atoms.atoms
    #         if isinstance(key, int):
    #             keys = [key]
    #         elif isinstance(key, tuple):
    #             keys = list(key)
    #         self.atoms[keys] = atoms
    #         self.symbols[keys] = atoms.symbols
    #         self.coords[keys] = atoms.coords
    #     else:
    #         self.atoms[key] = atoms
    #         self.symbols[key] = atoms.symbol
    #         self.coords[key] = atoms.coord

    def __setitem__(
        self, key: int | slice | list | tuple, atoms: Atom | Molecule
    ) -> None:
        if not isinstance(atoms, (Atom, Molecule)):
            raise TypeError(
                f"atoms must be an Atom or Molecule object, got {type(atoms).__name__}"
            )

        if isinstance(atoms, Molecule):
            if isinstance(key, int):
                keys = [key]
            elif isinstance(key, (list, tuple)):
                keys = list(key)
            elif isinstance(key, slice):
                keys = list(range(*key.indices(len(self.atoms))))
            else:
                raise TypeError(f"Unsupported key type: {type(key).__name__}")

            if len(keys) != len(atoms.atoms):
                raise ValueError(
                    f"Number of keys ({len(keys)}) does not match number of atoms in Molecule ({len(atoms.atoms)})"
                )

            for i, atom in zip(keys, atoms.atoms):
                self.atoms[i] = atom
                self.symbols[i] = atom.symbol
                self.coords[i] = atom.coord
        else:
            if isinstance(key, int):
                self.atoms[key] = atoms
                self.symbols[key] = atoms.symbol
                self.coords[key] = atoms.coord
            else:
                raise ValueError("Key must be an integer when assigning a single Atom.")

    def __iter__(self):
        return iter(self.atoms)

    def __str__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    def __repr__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    def sorted(
        self,
        key: callable = lambda atom: [
            ATOM_SYMBOLS[atom.symbol],
            atom.coord[0],
            atom.coord[1],
            atom.coord[2],
        ],
    ):
        indices = sorted(range(len(self.atoms)), key=lambda i: key(self.atoms[i]))
        return Molecule.from_atoms(
            [self.atoms[i] for i in indices], lattice_vecs=self.lattice_vecs
        )

    def sort(
        self,
        key: callable = lambda atom: [
            ATOM_SYMBOLS[atom.symbol],
            atom.coord[0],
            atom.coord[1],
            atom.coord[2],
        ],
    ):
        indices = sorted(range(len(self.atoms)), key=lambda i: key(self.atoms[i]))
        for i in indices:
            self.atoms[i].symbol = self.symbols[i]
            self.atoms[i].coord[:] = self.coords[i]
        self.symbols[:] = self.symbols[indices]
        self.coords[:] = self.coords[indices]

    def pop(self, index: int) -> Atom:
        return self.atoms.pop(index)

    def copy(self) -> Molecule:
        new_mol = Molecule(*copy.deepcopy(self.atoms))
        new_mol.lattice_vecs = copy.deepcopy(self.lattice_vecs)
        return new_mol

    def is_same_geom(
        self, other: Molecule, rel_tol: float = 1e-5, abs_tol: float = 0.0
    ) -> bool:
        sorted_self = self.sorted()
        sorted_other = other.sorted()
        return len(self) == len(other) and all(
            atom1.isclose(atom2, rel_tol=rel_tol, abs_tol=abs_tol)
            for atom1, atom2 in zip(
                sorted_self,
                sorted_other,
            )
        )

    @args_to_set
    def filtered_by_symbols(self, symbols: str | Iterable[str]) -> Molecule:
        return Molecule(*[atom for atom in self if atom.symbol in symbols])

    @classmethod
    @args_to_list
    def from_atoms(
        cls, atoms: Iterable[Atom], lattice_vecs: ArrayLike | None = None
    ) -> Molecule:
        """
        Create a new molecule from Atom objects.
        """
        if not all(isinstance(atom, Atom) for atom in atoms):
            raise TypeError("All elements must be Atom objects")
        return cls(*atoms, lattice_vecs=lattice_vecs)

    def add_atom(self, new_atom: Atom) -> None:
        """
        Add an atom to the molecule.
        makes sure self.symbols and self.coords share the same memory as the self.atoms
        """
        if not isinstance(new_atom, Atom):
            raise TypeError("new_atom must be an Atom object")

        oldN = len(self.atoms)

        if oldN == 0:
            self.coords = np.asarray([new_atom.coord], dtype=float)
            self.symbols = np.asarray([new_atom.symbol], dtype="<U2")
        else:
            new_coords = np.zeros((oldN + 1, 3), dtype=float)
            new_coords[:oldN, :] = self.coords
            new_coords[-1, :] = new_atom.coord

            new_symbols = np.zeros((oldN + 1,), dtype="<U2")
            new_symbols[:oldN] = self.symbols
            new_symbols[-1] = new_atom.symbol

            self.coords = new_coords
            self.symbols = new_symbols

        self.atoms.append(new_atom)

        # Update the atom coordinates to share the same memory
        for i, atom in enumerate(self.atoms):
            atom.symbol = self.symbols[i]
            atom.coord = self.coords[i]

    # def add_atoms_from(self, new_atoms: Iterable[Atom]) -> None:
    #     """
    #     Add atoms to the molecule.
    #     Accepts Iterable containing one or more Atom objects.
    #     """
    #     if any(not isinstance(atom, Atom) for atom in new_atoms):
    #         raise TypeError("All elements in atoms must be Atom objects")
    #
    #     oldN = len(self.atoms)
    #
    #     if oldN == 0:
    #         self.coords = np.asarray(
    #             [atom.coord for atom in new_atoms], dtype=np.float64
    #         )
    #         self.symbols = np.asarray([atom.symbol for atom in new_atoms], dtype="<U2")
    #     else:
    #         new_coords = np.zeros((oldN + len(new_atoms), 3), dtype=float)
    #         new_coords[:oldN, :] = self.coords
    #         new_coords[oldN:, :] = [atom.coord for atom in new_atoms]
    #
    #         new_symbols = np.zeros((oldN + len(new_atoms),), dtype="<U2")
    #         new_symbols[:oldN] = self.symbols
    #         new_symbols[oldN:] = [atom.symbol for atom in new_atoms]
    #
    #         self.coords = new_coords
    #         self.symbols = new_symbols
    #
    #     self.atoms.extend(new_atoms)
    #
    #     # Update the atom coordinates to share the same memory
    #     for i, atom in enumerate(self.atoms):
    #         atom.symbol = self.symbols[i]
    #         atom.coord = self.coords[i]

    def add_atoms_from(self, new_atoms: Iterable[Atom]) -> None:
        """
        Add atoms to the molecule and ensure memory sharing.
        """
        if any(not isinstance(atom, Atom) for atom in new_atoms):
            raise TypeError("All elements in atoms must be Atom objects")

        new_coords = np.array([atom.coord for atom in new_atoms], dtype=np.float64)
        new_symbols = np.array([atom.symbol for atom in new_atoms], dtype="<U2")

        if len(self.atoms) == 0:
            self.coords = new_coords
            self.symbols = new_symbols
        else:
            self.coords = np.vstack([self.coords, new_coords])
            self.symbols = np.hstack([self.symbols, new_symbols])

        self.atoms.extend(new_atoms)

        for i in range(len(self.atoms)):
            self.atoms[i].coord = self.coords[i]
            self.atoms[i].symbol = self.symbols[i]

    def _get_geom_hash(self):
        return hash(tuple((atom.symbol, atom.x, atom.y, atom.z) for atom in self))

    def get_formula(self) -> str:
        symbol_count = dict()
        for atom in self:
            symbol_count[atom.symbol] = symbol_count.get(atom.symbol, 0) + 1

        formula = "-".join(
            f"{symbol}{count}" if count > 1 else symbol
            for symbol, count in sorted(
                symbol_count.items(), key=lambda x: ATOM_SYMBOLS[x[0]]
            )
        )
        return formula

    def get_frac_coords(self, wrap=False) -> np.ndarray:
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to bound the molecule.")

        return cart2frac(
            cart_coords=self.coords, lattice_vecs=self.lattice_vecs, wrap=wrap
        )

    @cachedmethod(
        lambda self: self._cache,
        key=lambda self, tol=0.15: hashkey(self._get_geom_hash(), tol),
    )
    def get_bonds(
        self,
        tol: float = default_tol,
    ) -> list[dict[str, tuple[int, int], float]]:
        """
        Get the bonds in the molecule.
        caches the result for the same geometry and tolerance
        """

        bonds = list()
        num_atoms = len(self)
        for i in range(num_atoms):
            ai = self[i]
            for j in range(i + 1, num_atoms):
                aj = self[j]
                length = ai.get_bond_length(aj, tol)
                if length is not None:
                    bonds.append({"pair": (i, j), "length": length})

        return bonds

    def get_clusters(self, tol: float = default_tol) -> list[Molecule]:
        G = nx.Graph()
        G.add_nodes_from(range(len(self)))
        G.add_edges_from(bond["pair"] for bond in self.get_bonds(tol))
        copied_mol = self.copy()
        return [copied_mol[list(cluster)] for cluster in nx.connected_components(G)]

    def get_cycles(
        self, length_bound: int | None = None, tol: float = default_tol
    ) -> list[Molecule]:
        G = nx.Graph()
        G.add_edges_from(bond["pair"] for bond in self.get_bonds(tol))
        cycles = nx.simple_cycles(G, length_bound=length_bound)
        return [self[list(cycle)] for cycle in cycles]

    def get_connected_cluster(
        self, atom_idx: int, tol: float = default_tol
    ) -> Molecule:
        G = nx.Graph()
        G.add_edges_from(bond["pair"] for bond in self.get_bonds(tol))
        return self[list(nx.node_connected_component(G, atom_idx))]

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

    def merge(self, mols: Molecule | Iterable[Molecule]) -> Molecule:
        """
        Merge other Molecule objects into this Molecule object.
        """
        if isinstance(mols, Molecule):
            mols = [mols]
        elif isinstance(mols, Iterable):
            mols = list(mols)

        if not all(isinstance(mol, Molecule) for mol in mols):
            raise TypeError(
                "All elements must be Molecule objects\n"
                + f"{[(i, type(mol)) for i, mol in enumerate(mols) if not isinstance(mol, Molecule)]}"
            )

        for mol in mols:
            self.add_atoms_from(mol)

    def translate(
        self, trans_vec: Vec3 | ArrayLike, with_lattice_vecs: bool = False
    ) -> None:

        trans_vec = np.asarray(trans_vec, dtype=float)
        if trans_vec.shape != (3,):
            raise ValueError(
                f"trans_vec must be a 3 element array, got {trans_vec.shape}"
            )

        self.coords[:] = self.coords + trans_vec

        if with_lattice_vecs and self.lattice_vecs is not None:
            self.lattice_vecs[:] = self.lattice_vecs + trans_vec

    def mirror(self, sx: int, sy: int, sz: int) -> None:
        for atom in self.atoms:
            atom.mirror(sx, sy, sz)

    def mirror_by_plane(self, p1: Vec3, p2: Vec3, p3: Vec3) -> None:
        for atom in self.atoms:
            atom.mirror_by_plane(p1, p2, p3)

    def matmul(self, mat: ArrayLike, with_lattice_vecs: bool = True) -> None:
        mat = np.asarray(mat, dtype=float)
        if mat.shape != (3, 3):
            raise ValueError(f"mat must be a 3x3 matrix, got shape {mat.shape}")

        self.coords[:] = self.coords @ mat.T

        if with_lattice_vecs and self.lattice_vecs is not None:
            self.lattice_vecs = self.lattice_vecs @ mat

    def rotate_by_axis(
        self,
        axis_point1: VecLike,
        axis_point2: VecLike,
        deg: float | int,
        with_lattice_vecs: bool = True,
    ) -> None:
        """
        :param axis_point1: One point on the rotation axis
        :param axis_point2: Another point on the rotation axis
        :param angle_degrees: Rotation angle (degrees)
        """
        axis_point1 = np.asarray(axis_point1)
        axis_point2 = np.asarray(axis_point2)
        angle_radians = np.deg2rad(deg)

        axis_vector = axis_point2 - axis_point1
        axis_unit_vector = axis_vector / np.linalg.norm(axis_vector)
        rot = Rotation.from_rotvec(angle_radians * axis_unit_vector)
        point = self.coords - axis_point1
        rotated_point = rot.apply(point)

        self.coords[:] = rotated_point + axis_point1

        if with_lattice_vecs and self.lattice_vecs is not None:
            self.lattice_vecs[:] = (
                rot.apply(self.lattice_vecs - axis_point1) + axis_point1
            )

    def replicated(
        self, a_range: list[int, int], b_range: list[int, int], c_range: list[int, int]
    ) -> Molecule:
        """
        Replicate the molecule in the a, b, and c directions.
        """
        if len(self) == 0:
            raise ValueError("Molecule must have atoms to replicate, got 0 atoms")

        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to replicate the molecule.")

        if not all(isinstance(rep, list) for rep in (a_range, b_range, c_range)):
            raise TypeError(
                "Replication vectors must be of length 2 list\n"
                + f"{(a_range, b_range, c_range)=}"
            )

        if len(a_range) != 2 or len(b_range) != 2 or len(c_range) != 2:
            raise ValueError(
                "Replication vectors must be of length 2 list\n"
                + f"{(a_range, b_range, c_range)=}"
            )

        if (
            a_range[0] >= a_range[1]
            or b_range[0] >= b_range[1]
            or c_range[0] >= c_range[1]
        ):
            raise ValueError(
                "Replication vectors must be of the form [start, end], start < end\n"
                + f"{(a_range, b_range, c_range)=}"
            )

        new_mol = Molecule(lattice_vecs=self.lattice_vecs.copy())
        for i in range(a_range[0], a_range[1]):
            for j in range(b_range[0], b_range[1]):
                for k in range(c_range[0], c_range[1]):

                    trans_vec = (
                        i * new_mol.lattice_vecs[0]
                        + j * new_mol.lattice_vecs[1]
                        + k * new_mol.lattice_vecs[2]
                    )
                    mol_copied = self.copy()
                    mol_copied.translate(trans_vec)
                    new_mol.merge(mol_copied)

        return new_mol

    def is_inside_cell(self, atom: Atom) -> bool:
        """
        Check if the atom is inside the parallelepiped defined by the lattice vectors.

        let v = u1*a1 + u2*a2 + u3*a3, where a1, a2, a3 are lattice vectors
        mat = [a2×a3, a3×a1, a1×a2]
        to check if v is inside the cell, we need to check if 0 <= u1, u2, u3 < 1
        U = [u1, u2, u3] = 1/V * mat * v
        """
        if not isinstance(atom, Atom):
            raise TypeError(f"Expected Atom, got {type(atom).__name__}")

        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to bound the molecule.")

        U = atom.get_frac_coords(self.lattice_vecs, wrap=False)
        return np.all((0 <= U) & (U < 1))

    def remove_duplicates(self, tol: float = default_tol) -> None:
        """
        Remove duplicate atoms by combining close atoms with the same symbol
        into their centroid.
        Args:
            tol (float): Tolerance for distance to consider atoms as duplicates.
        """
        unique_coords = []
        unique_symbols = []
        unique_atoms = []

        symbols = np.unique(self.symbols)
        for symbol in symbols:
            indices = [i for i, sym in enumerate(self.symbols) if sym == symbol]
            sub_coords = self.coords[indices]
            sub_atoms = [self.atoms[i] for i in indices]

            # Use KDTree to cluster atoms within the tolerance
            tree = cKDTree(sub_coords)
            clusters = tree.query_ball_tree(tree, tol)

            processed = set()
            for cluster in clusters:
                if any(idx in processed for idx in cluster):
                    continue

                processed.update(cluster)

                cluster_coords = sub_coords[cluster]
                cluster_atoms = [sub_atoms[i] for i in cluster]
                # can calculate center of mass like this because every symbol is the same
                centroid = cluster_coords.mean(axis=0)

                unique_coords.append(centroid)
                unique_symbols.append(symbol)
                unique_atoms.append(Atom(symbol, *centroid))

        self.coords = np.array(unique_coords)
        self.symbols = np.array(unique_symbols)
        self.atoms = unique_atoms

    def wrap_to_cell(self, remove_duplicates: bool = False) -> None:
        """
        Wrap the molecule to the cell defined by the lattice vectors.
        """
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to bound the molecule.")

        self.coords[:] = frac2cart(self.get_frac_coords(wrap=True), self.lattice_vecs)

        if remove_duplicates:
            self.remove_duplicates()

    def replicated_from_xyz_str(
        self, xyz_str: str, wrap: bool = True, remove_duplicates: bool = False
    ) -> Molecule:
        """
        Replicate the molecule from an xyz string in cif format.
        args:
            xyz_str (str): xyz string of symmetry operation. e.g. 'x, y, z', '-y, -x+3/4, z+1/2',
            wrap (bool): wrap the molecule to the cell after replication
        """

        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to replicate the molecule.")

        rot_mat, trans_vec_frac = symmop_from_xyz_str(xyz_str)

        new_mol = self.copy()
        new_mol.coords[:] = new_mol.coords @ rot_mat.T
        new_mol.translate(frac2cart(trans_vec_frac, self.lattice_vecs))

        if wrap:
            new_mol.wrap_to_cell(remove_duplicates=remove_duplicates)

        return new_mol

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
        dist_mat = np.linalg.norm(
            self.coords[:, np.newaxis, :] - self.coords[np.newaxis, :, :], axis=-1
        )
        np.fill_diagonal(dist_mat, 1.0)
        charges = np.array([atom.atomic_number for atom in self])
        # sum up the upper triangle of the matrix to avoid double counting
        nuclrep = np.sum(np.triu(charges / dist_mat, k=1))
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
        dist_mat = np.linalg.norm(
            self.coords[:, np.newaxis, :] - self.coords[np.newaxis, :, :], axis=-1
        )
        np.fill_diagonal(dist_mat, 1.0)
        charges = []
        for i, atom in enumerate(self):
            if atom.charge is None:
                raise ValueError(f"Charge of atom {i} {atom} is None")
            charges.append(atom.charge)
        charges = np.array(charges)
        elec_energy = np.sum(np.triu(charges / dist_mat, k=1))
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

    def write_to_gaussian_input(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write("#p B3LYP\n")
            f.write("\n")
            f.write(f"{self.get_formula()}\n")
            f.write("\n")
            f.write("0 1\n")
            f.write(self.to_xyz() + "\n")
            if self.lattice_vecs is not None:
                for vec in self.lattice_vecs:
                    f.write(
                        f"{'Tv':2s} {vec[0]:19.12f} {vec[1]:19.12f} {vec[2]:19.12f}\n"
                    )
        print(f"File written to {filepath}")

    def write_to_gamess_input(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write("$DATA\n")
            f.write(f"{self.get_formula()}\n")
            f.write("C1\n")
            for i in range(len(self)):
                f.write(
                    f"{self.symbols[i]:2s}  {float(self[i].atomic_number):3d}  {self[i].x:19.12f} {self[i].y:19.12f} {self[i].z:19.12f}\n"
                )
        print(f"File written to {filepath}")

    def write_to_poscar(self, filepath: str, frac: bool = True) -> None:
        with open(filepath, "w") as f:
            f.write(f"{self.get_formula()}\n")
            f.write("1.0\n")
            f.write(
                f"{self.lattice_vecs[0][0]:19.12f} {self.lattice_vecs[0][1]:19.12f} {self.lattice_vecs[0][2]:19.12f}\n"
            )
            f.write(
                f"{self.lattice_vecs[1][0]:19.12f} {self.lattice_vecs[1][1]:19.12f} {self.lattice_vecs[1][2]:19.12f}\n"
            )
            f.write(
                f"{self.lattice_vecs[2][0]:19.12f} {self.lattice_vecs[2][1]:19.12f} {self.lattice_vecs[2][2]:19.12f}\n"
            )
            unique_symbols = np.unique(self.symbols)
            f.write(" ".join(unique_symbols) + "\n")
            symbol_count = {
                symbol: np.sum(self.symbols == symbol) for symbol in unique_symbols
            }
            f.write(
                " ".join(str(symbol_count[symbol]) for symbol in unique_symbols) + "\n"
            )
            if frac:
                f.write("Direct\n")
                coords = self.get_frac_coords(wrap=False)
            else:
                f.write("Cartesian\n")
                coords = self.coords
            for coord in coords:
                f.write(f"{coord[0]:19.12f} {coord[1]:19.12f} {coord[2]:19.12f}\n")
        print(f"File written to {filepath}")

    def write_to_cif(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write("data_molecule\n")
            f.write("########################\n")
            f.write("# generated by molgeom #\n")
            f.write("########################\n")
            f.wirte("\n")
            f.write("loop_\n")
            f.write("_space_group_symop_operation_xyz\n")
            f.write(" 'x, y, z'\n")
            f.write("\n")
            a, b, c, alpha, beta, gamma = lat_vecs_to_lat_params(self.lattice_vecs)
            f.write(f"_cell_length_a {a:.9f}\n")
            f.write(f"_cell_length_b {b:.9f}\n")
            f.write(f"_cell_length_c {c:.9f}\n")
            f.write(f"_cell_angle_alpha {alpha:.9f}\n")
            f.write(f"_cell_angle_beta {beta:.9f}\n")
            f.write(f"_cell_angle_gamma {gamma:.9f}\n")
            f.write("\n")
            f.write("loop_\n")
            f.write("_atom_site_label\n")
            f.write("_atom_site_type_symbol\n")
            f.write("_atom_site_fract_x\n")
            f.write("_atom_site_fract_y\n")
            f.write("_atom_site_fract_z\n")
            frac_coords = self.get_frac_coords(wrap=False)
            for i in range(len(self)):
                f.write(
                    f"{i+1:3d} {atom.symbols[i]:2s} {frac_coords[i][0]:19.12f} {frac_coords[i][1]:19.12f} {frac_coords[i][2]:19.12f}\n"
                )
        print(f"File written to {filepath}")
