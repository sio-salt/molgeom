import copy

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


def test_replicate():
    a1 = Atom("O", 0, 0, 0)
    a2 = Atom("H", 1, 0, 0)
    a3 = Atom("H", 0, 1, 0)
    mol = Molecule(a1, a2, a3)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mol.replicate([0, 1], [0, 1], [0, 1])

    assert len(mol.atoms) == 3
    assert mol == Molecule(a1, a2, a3)

    mol = Molecule(a1, a2, a3)
    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mol_copied = mol.copy()
    mol_copied.translate(
        -1 * (mol.lattice_vecs[0] + mol.lattice_vecs[1] + mol.lattice_vecs[2])
    )
    mol.replicate([-1, 0], [-1, 0], [-1, 0])

    assert len(mol.atoms) == 3
    assert mol == mol_copied
    assert mol.lattice_vecs is not None

    mol.lattice_vecs = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mol.replicate([-1, 2], [-1, 2], [-1, 2])
    assert len(mol.atoms) == 81


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
    assert mol.lattice_vecs == mol_replicated.lattice_vecs

    xyz_str_2 = "y, -x, z"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_2)
    tmp_mol = Molecule(Atom("O", 1, -1, 1), Atom("H", 0, -1, 0), Atom("H", 1, 0, 0))
    print(mol_replicated.atoms)
    print(tmp_mol.atoms)
    assert mol_replicated == tmp_mol
    assert mol_replicated.lattice_vecs == [
        Vec3(*[0, -2, 0]),
        Vec3(*[2, 0, 0]),
        Vec3(*[0, 0, 2]),
    ]

    xyz_str_3 = "-x, y + 1/2, -z + 1/2"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_3)
    tmp_mol = Molecule(Atom("O", -1, 2, 0), Atom("H", -1, 1, 1), Atom("H", 0, 2, 1))
    assert mol_replicated == tmp_mol

    xyz_str_4 = "-2y+1/2, 3x+1/2, z-y+1/2"
    mol_replicated = mol.replicated_from_xyz_str(xyz_str_4)
    tmp_mol = Molecule(Atom("O", -1, 4, 1), Atom("H", 1, 4, 1), Atom("H", -1, 1, 0))
