import os
from molgeom import read_file


def test_read_file():
    basepath = os.path.join(os.path.dirname(__file__), "files/")
    filenames = os.listdir(basepath)
    print(filenames)
    for filename in filenames:
        mole = read_file(basepath + filename)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0


def test_xyz_parser():
    basepath = os.path.join(os.path.dirname(__file__), "files/")
    filepnames = [files for files in os.listdir(basepath) if files.endswith(".xyz")]
    for filepname in filepnames:
        mole = read_file(basepath + filepname)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0
