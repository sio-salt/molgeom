from molgeom import Atom, Vec3


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

    assert type(a1) == type(a2) == type(a3) == Atom
    assert a1 == a2 == a3
    assert (
        a1.x == a2.x == a3.x == 2
        and a1.y == a2.y == a3.y == 1.1
        and a1.z == a2.z == a3.z == -0.2
    )


def test_copy():
    atom = Atom("C", 0, 1.1, -0.2)
    atom_copied = atom.copy()

    assert atom == atom_copied
    assert atom is not atom_copied
