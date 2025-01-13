from __future__ import annotations

import copy
from pathlib import Path
from collections.abc import Iterable

from cachetools import cachedmethod, LRUCache
from cachetools.keys import hashkey
import networkx as nx

from molgeom.data.consts import ANGST2BOHR_GAU16, ATOMIC_NUMBERS
from molgeom.utils.fancy_indexing_list import FancyIndexingList
from molgeom.utils.vec3 import Vec3, mat_type, vec_type
from molgeom.utils.mat3 import Mat3, is_mat_type
from molgeom.utils.decorators import args_to_list, args_to_set
from molgeom.utils.lattice_utils import cart2frac, frac2cart, lat_vecs_to_lat_params
from molgeom.utils.symmetry_utils import symmop_from_xyz_str
from molgeom.utils.html_utils import view_mol
from molgeom.atom import Atom


default_tol = 0.15


class Molecule:
    """
    A class to represent a molecule.
    """

    def __init__(self, *atoms: Atom, lattice_vecs: mat_type | None = None):
        self.atoms: FancyIndexingList[Atom] = FancyIndexingList()
        if atoms:
            self.add_atoms_from(atoms)
        self._lattice_vecs: list[Vec3] | None = None
        self.lattice_vecs = lattice_vecs
        self._cache = LRUCache(maxsize=10)

    @property
    def lattice_vecs(self) -> list[Vec3] | None:
        return self._lattice_vecs

    @lattice_vecs.setter
    def lattice_vecs(self, lattice_vecs: mat_type | None) -> None:
        if lattice_vecs is None:
            self._lattice_vecs = None
            return

        lat_vecs = Mat3(lattice_vecs)
        lat_det = Mat3.det(lat_vecs)
        if lat_det == 0:
            raise ValueError("lattice vectors must not be linearly dependent")
        if lat_det < 0:
            raise ValueError("lattice vectors must form a right-handed basis")

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

    def __setitem__(self, index: int, atom: Atom) -> None:
        self.atoms[index] = atom

    def __delitem__(self, idx: int) -> None:
        del self.atoms[idx]

    def __iter__(self):
        for atom in self.atoms:
            if not isinstance(atom, Atom):
                raise TypeError(
                    "Invalid element type: atom must be Atom object\n" + f"{type(atom)=}"
                )
            yield atom

    def __str__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    def __repr__(self):
        formula = self.get_formula()
        return f"Molecule({formula})"

    @classmethod
    def sorted(cls, atoms: Molecule | Iterable[Atom], key: callable = None) -> Molecule:
        """
        Create a new molecule by sorting the atoms of the input molecule(s).
        """
        if isinstance(atoms, Molecule):
            atoms = atoms.atoms
        if not all(isinstance(atoms, Atom) for atoms in atoms):
            raise TypeError("All elements must be Molecule objects")
        atoms.sort(key=key)
        mol = cls()
        mol.add_atoms_from(atoms)

        return mol

    def sort(self, key=None) -> None:
        """
        Sort the atoms of the molecule.
        """
        if key is None:
            self.atoms.sort(key=lambda atom: (ATOMIC_NUMBERS[atom.symbol], atom.x, atom.y, atom.z))
        else:
            self.atoms.sort(key=key)

    def pop(self, index: int) -> Atom:
        return self.atoms.pop(index)

    def copy(self) -> Molecule:
        return copy.deepcopy(self)

    def is_same_geom(self, other: Molecule, rel_tol: float = 1e-5, abs_tol: float = 0.0) -> bool:
        sorted_self = Molecule.sorted(self)
        sorted_other = Molecule.sorted(other)
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
    def from_atoms(cls, atoms: Iterable[Atom]) -> Molecule:
        """
        Create a new molecule from Atom objects.
        """
        if not all(isinstance(atom, Atom) for atom in atoms):
            raise TypeError("All elements must be Atom objects")
        return cls(*atoms)

    @classmethod
    def from_file(cls, filepath: str | Path) -> Molecule:
        """
        Create a new molecule from a file.
        auto-detects file format (*.xyz, *.com, *.gjf, *.inp, *.cif, *POSCAR*)
        """
        from molgeom.parsers.parser_selector import read_file

        return read_file(filepath)

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
            raise TypeError("Invalid element type: atoms must be Atom objects " + f"{not_atoms=}")
        self.atoms.extend(atoms)

    def _get_geom_hash(self):
        return hash(tuple((atom.symbol, atom.x, atom.y, atom.z) for atom in self))

    def get_symbols(self) -> list[str]:
        return tuple(atom.symbol for atom in self)

    def get_formula(self) -> str:
        symbol_count = dict()
        for atom in self:
            symbol_count[atom.symbol] = symbol_count.get(atom.symbol, 0) + 1

        formula = "-".join(
            f"{symbol}{count}" if count > 1 else symbol
            for symbol, count in sorted(symbol_count.items(), key=lambda x: ATOMIC_NUMBERS[x[0]])
        )
        return formula

    def get_cart_coords(self) -> list[Vec3]:
        return [[atom.x, atom.y, atom.z] for atom in self]

    def get_frac_coords(self, wrap=False) -> list[Vec3]:
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to bound the molecule.")

        return [cart2frac(atom.to_Vec3(), self.lattice_vecs, wrap=wrap) for atom in self]

    @cachedmethod(
        lambda self: self._cache,
        key=lambda self, tol=0.15: hashkey(self._get_geom_hash(), tol),
    )
    def get_bonds(
        self,
        tol: float = default_tol,
    ) -> list[dict[tuple[int, int], float]]:

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

    def get_connected_cluster(self, atom_idx: int, tol: float = default_tol) -> Molecule:
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

    def translate(
        self, trans_vec: Vec3 | list[float | int], with_lattice_vecs: bool = False
    ) -> None:
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

        if with_lattice_vecs and self.lattice_vecs is not None:
            for vec in self.lattice_vecs:
                vec.translate(trans_vec)

    def mirror(self, sx: int, sy: int, sz: int) -> None:
        for atom in self.atoms:
            atom.mirror(sx, sy, sz)

    def mirror_by_plane(self, p1: Vec3, p2: Vec3, p3: Vec3) -> None:
        for atom in self.atoms:
            atom.mirror_by_plane(p1, p2, p3)

    def matmul(self, mat: mat_type, with_lattice_vecs: bool = True) -> None:
        if not is_mat_type(mat):
            raise TypeError("mat must be a 3x3 matrix")

        for atom in self.atoms:
            atom.matmul(mat)

        if with_lattice_vecs and self.lattice_vecs is not None:
            self.lattice_vecs = self.lattice_vecs @ mat

    def rotate_by_axis(
        self,
        axis_point1: vec_type,
        axis_point2: vec_type,
        deg: float,
        with_lattice_vecs: bool = True,
    ) -> None:
        """
        :param axis_point1: One point on the rotation axis
        :param axis_point2: Another point on the rotation axis
        :param angle_degrees: Rotation angle (degrees)
        """

        for atom in self.atoms:
            atom.rotate_by_axis(axis_point1, axis_point2, deg)

        if with_lattice_vecs and self.lattice_vecs is not None:
            for vec in self.lattice_vecs:
                vec.rotate_by_axis(axis_point1, axis_point2, deg)

    def replicate(self, rep_a: list[int], rep_b: list[int], rep_c: list[int]) -> None:
        """
        Replicate the molecule in the a, b, and c directions.
        """
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to replicate the molecule.")

        if not all(isinstance(rep, list) for rep in (rep_a, rep_b, rep_c)):
            raise TypeError(
                "Replication vectors must be of length 2 list\n" + f"{(rep_a, rep_b, rep_c)=}"
            )

        if len(rep_a) != 2 or len(rep_b) != 2 or len(rep_c) != 2:
            raise ValueError(
                "Replication vectors must be of length 2 list\n" + f"{(rep_a, rep_b, rep_c)=}"
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

        self.lattice_vecs = tmp_mol.lattice_vecs

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

        U = atom.get_frac_coords(self.lattice_vecs)
        return all(0 <= u < 1 for u in U)

    # def remove_duplicates(self, tol: float = default_tol) -> None:
    #     """
    #     Remove duplicate atoms (close atoms) from the molecule.
    #     finds too close clusters of same elements
    #     and combines them into one atom at the center of mass
    #     """
    #     clusters = self.get_clusters(tol)
    #     for cluster in clusters:
    #         if len(cluster) > 1:
    #             com = cluster.center_of_mass()
    #             cluster[:] = [Atom("C", com.x, com.y, com.z)]

    def remove_duplicates(self, tol: float = default_tol) -> None:
        """
        Remove duplicate atoms by combining close atoms with the same symbol
        into their centroid.
        Args:
            tol (float): Tolerance for distance to consider atoms as duplicates.
        """
        unique_atoms = FancyIndexingList()

        symbols = set(a.symbol for a in self.atoms)
        for symbol in symbols:
            indices = [i for i, atom in enumerate(self) if atom.symbol == symbol]
            sub_atoms = self[indices]

            processed = set()
            for i, atom in enumerate(sub_atoms):
                if i in processed:
                    continue

                cluster = [i]
                for j, other_atom in enumerate(sub_atoms):
                    if j != i and atom.distance_to(other_atom) < tol:
                        cluster.append(j)

                processed.update(cluster)

                cluster_atoms = sub_atoms[cluster]
                centroid = cluster_atoms.center_of_mass()

                unique_atoms.append(Atom(symbol, *centroid))

        self.atoms = unique_atoms

    def wrap_to_cell(self) -> None:
        """
        Wrap the molecule to the cell defined by the lattice vectors.
        """
        if self.lattice_vecs is None:
            raise ValueError("Lattice vectors must be set to bound the molecule.")

        frac_coords = self.get_frac_coords(wrap=True)
        for atom, frac_coord in zip(self, frac_coords):
            cart_coords = frac2cart(frac_coord, self.lattice_vecs)
            atom.x, atom.y, atom.z = cart_coords

    def replicated_from_xyz_str(self, xyz_str: str, wrap: bool = True) -> Molecule:
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
        new_mol.matmul(rot_mat, with_lattice_vecs=False)
        new_mol.translate(frac2cart(trans_vec_frac, self.lattice_vecs))

        if wrap:
            new_mol.wrap_to_cell()

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
        nuclrep = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                dist_angst = self.atoms[i].distance_to(self.atoms[j])
                dist_bohr = dist_angst * ANGST2BOHR_GAU16
                nuclrep += self.atoms[i].atomic_number * self.atoms[j].atomic_number / dist_bohr
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

    @staticmethod
    def write_to_xyzs(
        mols: list[Molecule],
        filepath: str | Path,
    ) -> None:
        with open(filepath, "w") as f:
            for mol in mols:
                f.write(f"{len(mol)}\n")
                f.write(mol.get_formula() + "\n")
                f.write(mol.to_xyz() + "\n")
        print(f"File written to {filepath}")

    def write_to_gaussian_input(
        self, filepath: str, head: str | None = None, tail: str | None = None
    ) -> None:
        with open(filepath, "w") as f:
            if head is not None:
                f.write(head)
            else:
                f.write("#p B3LYP\n")
                f.write("\n")
                f.write(f"{self.get_formula()}\n")
                f.write("\n")
                f.write("0 1\n")
            f.write(self.to_xyz() + "\n")
            if self.lattice_vecs is not None:
                for vec in self.lattice_vecs:
                    f.write(f"{'Tv':2s} {vec[0]:19.12f} {vec[1]:19.12f} {vec[2]:19.12f}\n")
            if tail is not None:
                f.write(tail)
        print(f"File written to {filepath}")

    def write_to_gamess_input(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            f.write("$DATA\n")
            f.write(f"{self.get_formula()}\n")
            f.write("C1\n")
            for i in range(len(self)):
                f.write(
                    f"{self[i].symbol:2s}  {float(self[i].atomic_number):1f}  {self[i].x:19.12f} {self[i].y:19.12f} {self[i].z:19.12f}\n"
                )
        print(f"File written to {filepath}")

    def write_to_poscar(self, filepath: str, frac: bool = True, wrap: bool = False) -> None:
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
            # remove duplicates and keep the order
            unique_symbols = list(dict.fromkeys([atom.symbol for atom in self]))
            f.write(" ".join(unique_symbols) + "\n")
            symbol_count = dict()
            for symbol in unique_symbols:
                symbol_count[symbol] = sum(atom.symbol == symbol for atom in self)
            f.write(" ".join(str(symbol_count[symbol]) for symbol in unique_symbols) + "\n")
            if frac:
                f.write("Direct\n")
                coords = self.get_frac_coords(wrap=wrap)
            else:
                f.write("Cartesian\n")
                coords = self.coords
            for coord in coords:
                f.write(f"{coord[0]:19.12f} {coord[1]:19.12f} {coord[2]:19.12f}\n")
        print(f"File written to {filepath}")

    def write_to_cif(self, filepath: str, wrap: bool = False) -> None:
        with open(filepath, "w") as f:
            f.write("data_molecule\n")
            f.write("########################\n")
            f.write("# generated by molgeom #\n")
            f.write("########################\n")
            f.write("\n")
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
            frac_coords = self.get_frac_coords(wrap=wrap)
            for i in range(len(self)):
                f.write(
                    f"{i+1:3d} {self[i].symbol:2s} {frac_coords[i][0]:19.12f} {frac_coords[i][1]:19.12f} {frac_coords[i][2]:19.12f}\n"
                )
        print(f"File written to {filepath}")

    def write_to_mol(self, filepath: str) -> None:
        """Write molecule to MOL V2000 format file.

        Creates a MOL file following V2000 format specification.
        Includes atomic coordinates, bonds, and charges.
        """
        bonds = self.get_bonds()

        with open(filepath, "w") as f:
            # Header block
            f.write(f"{self.get_formula()}\n")  # Molecule name
            f.write("  Generated by molgeom\n")  # Program info
            f.write("\n")  # Comment line

            # Counts line with actual bond count
            f.write(f"{len(self):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000\n")

            # Atom block
            for atom in self:
                charge = atom.charge if atom.charge is not None else 0
                f.write(
                    f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} {atom.symbol:<3s}"
                    f" 0  0  0  0  0  0  0  0  0{charge:3d}  0  0  0\n"
                )

            # Bond block - assuming single bonds for now
            for bond in bonds:
                i, j = bond["pair"]
                # MOL format uses 1-based indexing
                f.write(f"{i+1:3d}{j+1:3d}  1  0  0  0  0\n")

            # End marker
            f.write("M  END\n")

        print(f"File written to {filepath}")

    @staticmethod
    def write_to_sdf(
        mols: list[Molecule],
        filepath: str | Path,
        properties: dict[str, dict[str, str]] | None = None,
    ) -> None:
        """Write multiple Molecule objects to SDF format file.

        Args:
            mols: list of Molecule objects to write
            filepath: Output file path
            properties: Optional dictionary of molecular properties
                       Format: {mol_index: {property_name: property_value}}
        """
        with open(filepath, "w") as f:
            for i, mol in enumerate(mols):
                bonds = mol.get_bonds()

                # Write MOL format content
                f.write(f"{mol.get_formula()}\n")
                f.write("  Generated by molgeom\n")
                f.write("\n")
                f.write(f"{len(mol):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000\n")

                # Atom block
                for atom in mol:
                    charge = atom.charge if atom.charge is not None else 0
                    f.write(
                        f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} {atom.symbol:<3s}"
                        f" 0  0  0  0  0  0  0  0  0{charge:3d}  0  0  0\n"
                    )

                # Bond block
                for bond in bonds:
                    i, j = bond["pair"]
                    digits = len(str(len(mol)))
                    f.write(f"{i+1:>{digits+2}d}{j+1:>{digits+2}d}  1  0  0  0  0\n")

                f.write("M  END\n")

                # Write properties if provided
                if properties and i in properties:
                    for name, value in properties[i].items():
                        f.write(f"> <{name}>\n")
                        f.write(f"{value}\n")
                        f.write("\n")

                # Add molecule separator
                f.write("$$$$\n")

        print(f"File written to {filepath}")

    def show(self, cleanup: bool = True, prefer_notebook: bool = True) -> None:
        """
        View molecular geometry using 3Dmol.js in a browser.
        args:
            cleanup: bool
                If True, removes temporary files after viewing (default: True)
            prefer_notebook: bool
                If True, opens in Jupyter notebook if available (default: True)
        """
        xyz_data = f"{len(self)}\n{str(self)}\n{self.to_xyz()}"
        view_mol(xyz_mol_data=xyz_data, cleanup=cleanup, prefer_notebook=prefer_notebook)

    @staticmethod
    def show(mols: list[Molecule], cleanup: bool = True, prefer_notebook: bool = True) -> None:
        """
        View multiple molecular geometries using 3Dmol.js in a browser.
        args:
            mols: list[Molecule]
                List of Molecule objects to view
            cleanup: bool
                If True, removes temporary files after viewing (default: True)
            prefer_notebook: bool
                If True, opens in Jupyter notebook if available (default: True)
        """
        xyz_data = "\n".join([f"{len(mol)}\n{str(mol)}\n{mol.to_xyz()}" for mol in mols])
        view_mol(xyz_mol_data=xyz_data, cleanup=cleanup, prefer_notebook=prefer_notebook)
