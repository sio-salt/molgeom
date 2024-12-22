import copy

import numpy as np

from molgeom import Atom, Molecule, Vec3


def test_FancyIndexing_Molecule():
    a1 = Atom("N", 1.0, 1.0, 1.0)
    a2 = Atom("H", 1.0, 0.0, 0.0)
    a3 = Atom("H", 0.0, 1.0, 0.0)
    a4 = Atom("H", 0.0, 0.0, 1.0)
    mol = Molecule(a1, a2, a3, a4)

    assert mol[0] == a1
    assert mol[-1] == a4
    assert mol[[0, 1]] == Molecule(a1, a2)
    assert mol[(0, 1)] == Molecule(a1, a2)
    assert mol[0, 1] == Molecule(a1, a2)
    assert mol[0, -1] == Molecule(a1, a4)
    assert mol[[0, 1, 2]] == Molecule(a1, a2, a3)
    assert mol[[0, -1, 1, -2]] == Molecule(a1, a4, a2, a3)
    assert mol[1:] == Molecule(a2, a3, a4)
    assert mol[:3] == Molecule(a1, a2, a3)
    assert mol[-2:] == Molecule(a3, a4)
    assert mol[:-2] == Molecule(a1, a2)
    assert mol[1:3] == Molecule(a2, a3)
    assert mol[:] == Molecule(a1, a2, a3, a4)
    assert mol[1:3:2] == Molecule(a2)
    assert mol[1:4:2] == Molecule(a2, a4)
    assert mol[1::2] == Molecule(a2, a4)
    assert mol[:3:2] == Molecule(a1, a3)
    assert mol[::2] == Molecule(a1, a3)
    assert mol[::-1] == Molecule(a4, a3, a2, a1)
    mol[1] = Atom("C", 1.0, 1.0, 1.0)
    assert mol[1] == Atom("C", 1.0, 1.0, 1.0)
    mol[[1, 2]] = Molecule(Atom("C", 1.0, 1.0, 1.0), Atom("C", 1.0, 1.0, 1.0))
    assert mol[[1, 2]] == Molecule(Atom("C", 1.0, 1.0, 1.0), Atom("C", 1.0, 1.0, 1.0))
    mol[1, -1, 2] = Molecule(
        Atom("C", 1.0, 0.0, 0.0),
        Atom("C", 0.0, 1.0, 0.0),
        Atom("C", 0.0, 0.0, 1.0),
    )
    assert mol[1, -1, 2] == Molecule(
        Atom("C", 1.0, 0.0, 0.0),
        Atom("C", 0.0, 1.0, 0.0),
        Atom("C", 0.0, 0.0, 1.0),
    )
    assert mol[1:3] == Molecule(Atom("C", 1.0, 0.0, 0.0), Atom("C", 0.0, 0.0, 1.0))
    assert mol[0:4:2] == Molecule(Atom("N", 1.0, 1.0, 1.0), Atom("C", 0.0, 0.0, 1.0))
    mol = Molecule(a1, a2, a3, a4)
    mol[[0, 1, -1]] = Molecule(
        Atom("C", 1.0, 0.0, 0.0),
        Atom("C", 0.0, 1.0, 0.0),
        Atom("C", 0.0, 0.0, 1.0),
    )
    assert mol[[0, 1, -1]] == Molecule(
        Atom("C", 1.0, 0.0, 0.0),
        Atom("C", 0.0, 1.0, 0.0),
        Atom("C", 0.0, 0.0, 1.0),
    )
    assert mol == Molecule(
        Atom("C", 1.0, 0.0, 0.0),
        Atom("C", 0.0, 1.0, 0.0),
        Atom("H", 0.0, 1.0, 0.0),
        Atom("C", 0.0, 0.0, 1.0),
    )


def test_slice_id_check_Molecule():
    a1 = Atom("N", 1.0, 1.0, 1.0)
    a2 = Atom("H", 1.0, 0.0, 0.0)
    a3 = Atom("H", 0.0, 1.0, 0.0)
    a4 = Atom("H", 0.0, 0.0, 1.0)
    mol = Molecule(a1, a2, a3, a4)

    assert mol[:] == mol and mol[:] is not mol
    assert mol[1:3] == Molecule(a2, a3) and mol[1:3] is not Molecule(a2, a3)


def test_memory_share():
    a1 = Atom("N", 1.0, 1.0, 1.0)
    a2 = Atom("H", 1.0, 0.0, 0.0)
    a3 = Atom("H", 0.0, 1.0, 0.0)
    a4 = Atom("H", 0.0, 0.0, 1.0)
    mol = Molecule()
    mol.add_atoms_from([a1, a2, a3, a4])
    na1, na2, na3, na4 = mol[0], mol[1], mol[2], mol[3]

    assert mol.coords is mol._coords
    assert mol.symbols is mol._symbols

    mol.coords[2, 0] = 30.0
    assert mol[2].coord[0] == 30.0
    assert mol.coords[2, 0] == 30.0

    # Modify atom coordinate and check propagation
    a1.coord[0] = 10.0
    assert mol[0].coord[0] == 10.0
    assert mol.atoms[0].coord[0] == 10.0
    assert mol.coords[0][0] == 10.0
    assert mol.coords[0, 0] == 10.0

    a1.x = -100.0
    assert mol[0].coord[0] == -100.0 and mol.coords[0, 0] == -100.0

    mol[1].coord[0] = 20.0
    assert na2.coord[0] == 20.0 and mol[1].coord[0] == 20.0 and mol.coords[1, 0] == 20.0

    mol.coords[2, 0] = 30.0
    assert na3.coord[0] == 30.0 and mol[2].coord[0] == 30.0 and mol.coords[2, 0] == 30.0

    mol.coords[:] = mol.coords * 2
    assert mol[0].coord[0] == -200.0
    assert mol.atoms[0].coord[0] == -200.0
    assert mol[1].coord[0] == 40.0
    assert mol[2].coord[0] == 60.0
    assert mol[3].coord[0] == 0.0

    mol_coords = mol.coords.copy()

    mol.translate([1, 1, 1])
    assert np.array_equal(mol.coords, np.array([atom.coord.copy() for atom in mol]))

    angle = np.radians(30)
    mol.matmul(
        np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )
    )
    assert mol.coords[0, 0] == mol.atoms[0].coord[0] == mol.atoms[0].x
    assert np.array_equal(mol.coords, np.array([atom.coord.copy() for atom in mol]))

    mol.coords[:] = mol_coords
    assert np.array_equal(mol.coords, np.array([atom.coord.copy() for atom in mol]))

    mol.rotate_by_axis(Vec3(1, 1, 1), Vec3(2, 1, 1), deg=90.0)


def test_add_atom():
    mol1 = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    mol1.add_atom(a1)
    mol1.add_atoms_from([a2, a3])
    mol2 = Molecule(a1, a2, a3)

    assert len(mol1.atoms) == 3
    assert mol1.atoms[0] == a1
    assert mol1.atoms[1] == a2
    assert mol1.atoms[2] == a3
    assert mol1 == mol2


def test_copy():
    mol = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    mol.add_atoms_from([a1, a2, a3])
    mol_copied_from_method = mol.copy()
    mol_copied_from_copy = copy.deepcopy(mol)
    mol2 = mol

    assert mol == mol2
    assert mol is mol2
    assert mol == mol_copied_from_method
    assert mol is not mol_copied_from_method
    assert mol == mol_copied_from_copy
    assert mol is not mol_copied_from_copy


def test_is_same_geom():
    mol1 = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    mol1.add_atoms_from([a1, a2, a3])

    mol2 = Molecule()
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    a1 = Atom("C", 0, 1.1, 0.2)
    mol2.add_atoms_from([a1, a2, a3])

    assert mol1.is_same_geom(mol2)

    mol2 = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.5)
    mol2.add_atoms_from([a1, a2, a3])

    assert not mol1.is_same_geom(mol2)


def test_translate():
    mol = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    mol.add_atoms_from([a1, a2, a3])

    mol.translate([1, 1, 1])
    assert np.allclose(
        mol.coords, np.array([[1, 2.1, 1.2], [1, 2.2, 1.3], [1, 2.3, 1.4]])
    )


def test_rotate_by_axis():
    mol = Molecule()
    a1 = Atom("C", 0, 0, 1)
    mol.add_atoms_from([a1])

    mol.rotate_by_axis(Vec3(1, 0, 0), Vec3(1, 0, 1), deg=180.0)
    atom_coords = np.array([atom.coord.copy() for atom in mol])
    assert np.allclose(mol.coords, [[2, 0, 1]])
    assert np.allclose(mol.coords, atom_coords)


def test_matmul():
    mol = Molecule()
    a1 = Atom("C", 0, 1.1, 0.2)
    a2 = Atom("H", 0, 1.2, 0.3)
    a3 = Atom("H", 0, 1.3, 0.4)
    mol.add_atoms_from([a1, a2, a3])

    angle = np.radians(30)
    tmp = mol.copy()
    mat = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    mol.matmul(mat)
    assert np.allclose(mol.coords, tmp.coords @ mat.T)
    mol_atoms_coords = np.array([atom.coord.copy() for atom in mol])
    assert np.allclose(mol.coords, mol_atoms_coords)


def test_replicated():
    a1 = Atom("O", 0, 0, 0)
    a2 = Atom("H", 1, 0, 0)
    a3 = Atom("H", 0, 1, 0)
    mol = Molecule(a1, a2, a3)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    rep_mol = mol.replicated([0, 1], [0, 1], [0, 1])

    assert len(rep_mol.atoms) == 3
    assert rep_mol == Molecule(a1, a2, a3)

    mol = Molecule(a1, a2, a3)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mol_copied = mol.copy()
    mol_copied.translate(
        -1 * (mol.lattice_vecs[0] + mol.lattice_vecs[1] + mol.lattice_vecs[2])
    )
    rep_mol = mol.replicated([-1, 0], [-1, 0], [-1, 0])

    assert len(rep_mol.atoms) == 3
    assert rep_mol == mol_copied
    assert rep_mol.lattice_vecs is not None

    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    rep_mol = mol.replicated([-1, 2], [-1, 2], [-1, 2])
    assert len(rep_mol.atoms) == 81


def test_is_inside_cell():
    a1 = Atom("N", 3, 3, 3)
    a2 = Atom("H", 1, 0, 0)
    a3 = Atom("H", 0, 1, 0)
    a4 = Atom("H", 0, 0, 1)
    mol = Molecule(a1, a2, a3, a4)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    for atom in mol[1:]:
        assert mol.is_inside_cell(atom)
    assert not mol.is_inside_cell(mol[0])
    assert not mol.is_inside_cell(a1)

    mol = Molecule()
    mol.lattice_vecs = np.array([[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.3, 0.2, 1.0]])
    inside_vec = (
        0.1 * mol.lattice_vecs[0]
        + 0.2 * mol.lattice_vecs[1]
        + 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_1 = (
        1.1 * mol.lattice_vecs[0]
        + 0.2 * mol.lattice_vecs[1]
        + 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_2 = (
        0.1 * mol.lattice_vecs[0]
        + 1.2 * mol.lattice_vecs[1]
        + 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_3 = (
        0.1 * mol.lattice_vecs[0]
        + 0.2 * mol.lattice_vecs[1]
        + 1.3 * mol.lattice_vecs[2]
    )
    outside_vec_4 = (
        -1.5 * mol.lattice_vecs[0]
        + 0.2 * mol.lattice_vecs[1]
        + 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_5 = (
        0.1 * mol.lattice_vecs[0]
        - 0.2 * mol.lattice_vecs[1]
        + 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_6 = (
        0.1 * mol.lattice_vecs[0]
        + 0.2 * mol.lattice_vecs[1]
        - 0.3 * mol.lattice_vecs[2]
    )
    outside_vec_7 = (
        1.1 * mol.lattice_vecs[0]
        + 1.2 * mol.lattice_vecs[1]
        + 1.3 * mol.lattice_vecs[2]
    )
    outside_vec_8 = (
        -1.1 * mol.lattice_vecs[0]
        - 1.2 * mol.lattice_vecs[1]
        - 1.3 * mol.lattice_vecs[2]
    )
    assert mol.is_inside_cell(Atom.from_vec("H", inside_vec))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_1))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_2))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_3))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_4))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_5))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_6))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_7))
    assert not mol.is_inside_cell(Atom.from_vec("H", outside_vec_8))


def test_bound_to_cell():
    a1 = Atom("N", 3, 3, 3)
    a2 = Atom("H", 1, 0, 0)
    a3 = Atom("H", 0, 1, 0)
    a4 = Atom("H", 0, 0, 1)
    mol = Molecule(a1, a2, a3, a4)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

    mol.wrap_to_cell()
    for atom in mol:
        assert mol.is_inside_cell(atom)
    assert mol[0] == Atom("N", 1, 1, 1)


def test_replicated_from_xyz_str():
    a1 = Atom("O", 1, 1, 1)
    a2 = Atom("H", 1, 0, 0)
    a3 = Atom("H", 0, 1, 0)
    mol = Molecule(a1, a2, a3)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

    # e.g.   ‘x, y, z’, ‘-x, -y, z’, '-x, y + 1/2, -z + 1/2', ‘-2y+1/2, 3x+1/2, z-y+1/2’,
    xyz_str_1 = "x, y, z"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_1)
    assert mol == mol_replicated
    assert np.array_equal(mol.lattice_vecs, mol_replicated.lattice_vecs)

    mol = Molecule(Atom("O", 0.6, 1.1, 1.8))
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    xyz_str_2 = "y, -x, z"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_2)
    tmp_mol = Molecule(Atom("O", 1.1, 1.4, 1.8))
    assert mol_replicated == tmp_mol

    xyz_str_3 = "-x, y + 1/2, -z + 1/2"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_3)
    tmp_mol = Molecule(
        Atom("O", 1.4, 0.1, 1.2),
    )
    assert mol_replicated.is_same_geom(tmp_mol)

    xyz_str_4 = "-2y+1/2, 3x+1/2, z-y+1/2"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_4)
    tmp_mol = Molecule(Atom("O", 0.8, 0.8, 1.7))
    assert mol_replicated.is_same_geom(tmp_mol)
