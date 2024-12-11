from molgeom import read_file


def test_read_file():
    basepath = "./files/"
    filepnames = [
        "./H2O.xyz",
        "./1hvc.xyz",
        "./C14N2H8S4.xyz",
        "./Cr2C4H4O8.xyz",
        "./H2O_with_Tv.com",
        "./graphite_3layers.gjf",
        "./POSCAR_H2O",
    ]
    for filepname in filepnames:
        mole = read_file(basepath + filepname)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0


def test_xyz_parser():
    basepath = "./files/"
    filepnames = [
        "./H2O.xyz",
        "./1hvc.xyz",
        "./C14N2H8S4.xyz",
        "./Cr2C4H4O8.xyz",
    ]
    for filepname in filepnames:
        mole = read_file(basepath + filepname)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0
