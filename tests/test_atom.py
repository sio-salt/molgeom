import pytest

from molgeom import Atom, Vec3, Mat3


def test_instanciation():
    atom = Atom("C", 0, 1.1, -0.2)
    assert atom.symbol == "C"
    assert atom.x == 0
    assert atom.y == 1.1
    assert atom.z == -0.2


def test_from_vec():
    list_vec = [2, 1.1, -0.2]
    tuple_vec = (2, 1.1, -0.2)
    vec3 = Vec3(2, 1.1, -0.2)
    a1 = Atom.from_vec("C", list_vec)
    a2 = Atom.from_vec("C", tuple_vec)
    a3 = Atom.from_vec("C", vec3)
    assert isinstance(a1, Atom) and isinstance(a2, Atom) and isinstance(a3, Atom)
    assert a1 == a2 == a3
    assert (
        a1.x == a2.x == a3.x == 2
        and a1.y == a2.y == a3.y == 1.1
        and a1.z == a2.z == a3.z == -0.2
    )


def test_to_Vec3():
    atom = Atom("C", 0, 1.1, -0.2)
    vec = atom.to_Vec3()
    assert vec == Vec3(0, 1.1, -0.2)


def test_copy():
    atom = Atom("C", 0, 1.1, -0.2)
    atom_copied = atom.copy()
    assert atom == atom_copied
    assert atom is not atom_copied


def test_get_frac_coords():
    atom = Atom("C", 0.3, 1.1, -0.2)
    lattice_vec = Mat3([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    frac_coords = atom.get_frac_coords(lattice_vec, wrap=False)
    assert frac_coords.isclose(Vec3(0.3, 1.1, -0.2))
    frac_coords = atom.get_frac_coords(lattice_vec, wrap=True)
    assert frac_coords.isclose(Vec3(0.3, 0.1, 0.8))

    lattice_vec = Mat3([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    frac_coords = atom.get_frac_coords(lattice_vec, wrap=False)
    assert frac_coords.isclose(Vec3(0.15, 0.55, -0.1))
    frac_coords = atom.get_frac_coords(lattice_vec, wrap=True)
    assert frac_coords.isclose(Vec3(0.15, 0.55, 0.9))

    lattice_vec = Mat3([[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    frac_coords = atom.get_frac_coords(lattice_vec, wrap=True)
    assert frac_coords.isclose(Vec3(0.3, 0.55, 0.9))


def test_invalid_symbol():
    with pytest.raises(ValueError):
        Atom("no", 0, 0, 0)


def test_to_dict():
    atom = Atom("C", 0, 1.1, -0.2)
    atom_dict = atom.to_dict()
    assert atom_dict["symbol"] == "C"
    assert atom_dict["x"] == 0
    assert atom_dict["y"] == 1.1
    assert atom_dict["z"] == -0.2
    assert atom_dict["mass"] == atom.mass
    assert atom_dict["atomic_number"] == atom.atomic_number


def test_to_xyz():
    atom = Atom("C", 0, 1.1, -0.2)
    xyz = atom.to_xyz()
    assert xyz == "C       0.000000000000      1.100000000000     -0.200000000000"


def test_equality():
    atom1 = Atom("C", 0, 1.1, -0.2)
    atom2 = Atom("C", 0, 1.1, -0.2)
    atom3 = Atom("H", 0, 1.1, -0.2)
    assert atom1 == atom2
    assert atom1 != atom3


def test_get_bond_length():
    atom1 = Atom("C", 0, 0, 0)
    atom2 = Atom("C", 1.54, 0, 0)  # Typical C-C bond length
    bond_length = atom1.get_bond_length(atom2)
    assert bond_length is not None
    assert bond_length == pytest.approx(1.54, 0.01)
