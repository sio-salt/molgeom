import gzip
import bz2
import shutil
from pathlib import Path
from molgeom import read_file, Molecule
from molgeom.parsers import (
    extract_head_tail_from_gau_inp,
    extract_head_tail_from_gms_inp,
)


def test_read_file():
    basepath = Path(__file__).parent / "files"
    filenames = [f for f in basepath.iterdir() if f.is_file()]

    for filename in filenames:
        if filename.suffix in [".gz", ".bz2", ".xz", ".lzma"]:
            continue

        mole = read_file(filename)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0

        mole_2 = Molecule.from_file(filename)
        assert mole == mole_2

        file_suffix = filename.suffix

        gz_filename = filename.with_suffix(f"{file_suffix}.gz")
        with open(filename, "rb") as f_in:
            with gzip.open(gz_filename, "wb") as f_gz:
                shutil.copyfileobj(f_in, f_gz)

        bz2_filename = filename.with_suffix(f"{file_suffix}.bz2")
        with open(filename, "rb") as f_in:
            with bz2.open(bz2_filename, "wb") as f_bz2:
                shutil.copyfileobj(f_in, f_bz2)

        try:
            mole_gz = read_file(gz_filename)
            assert mole == mole_gz
        finally:
            if gz_filename.exists():
                gz_filename.unlink()

        try:
            mole_bz2 = read_file(bz2_filename)
            assert mole == mole_bz2
        finally:
            if bz2_filename.exists():
                bz2_filename.unlink()


def test_xyz_parser():
    basepath = Path(__file__).parent / "files"
    filepnames = [file for file in basepath.iterdir() if file.suffix == ".xyz"]
    for filepname in filepnames:
        mole = read_file(basepath / filepname)
        assert mole is not None
        assert mole.atoms is not None
        assert len(mole.atoms) > 0


def test_extract_head_tail_from_gau_inp():
    gau_file = Path(__file__).parent / "files" / "CH2_opt.com"
    head, tail = extract_head_tail_from_gau_inp(gau_file)
    file_head = """#p UHF/Gen Scf(MaxCycle=250)
opt=(z-matrix, maxcycle=30)
GFOldPrint pop=full Guess=Huckel

comment line

0 3"""
    file_tail = """
H     0
S    3   1.00
      0.1873113696D+02       0.3349460434D-01
      0.2825394365D+01       0.2347269535D+00
      0.6401216923D+00       0.8137573261D+00
S    1   1.00
      0.1612777588D+00       1.0000000
****
C     0
S    6   1.00
      0.3047524880D+04       0.1834737132D-02
      0.4573695180D+03       0.1403732281D-01
      0.1039486850D+03       0.6884262226D-01
      0.2921015530D+02       0.2321844432D+00
      0.9286662960D+01       0.4679413484D+00
      0.3163926960D+01       0.3623119853D+00
SP   3   1.00
      0.7868272350D+01      -0.1193324198D+00       0.6899906659D-01
      0.1881288540D+01      -0.1608541517D+00       0.3164239610D+00
      0.5442492580D+00       0.1143456438D+01       0.7443082909D+00
SP   1   1.00
      0.1687144782D+00       0.1000000000D+01       0.1000000000D+01
D    1   1.00
      0.8000000000D+00       1.0000000
****
"""
    assert head == file_head
    assert tail == file_tail


def test_extract_head_tail_from_gms_inp():
    gms_file = Path(__file__).parent / "files" / "I2_mp2.inp"
    head, tail = extract_head_tail_from_gms_inp(gms_file)
    file_head = """ $SYSTEM MWORDS= 2602 TIMLIM=2000000 $END
 $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG=0 MULT=1
         ISPHER=0 MAXIT=80 EXETYP=RUN PP=READ MPLEVL=2  $END
 $BASIS  BASNAM(1)=I                                    $END
 $SCF    DIRSCF=.F. DIIS=.T. SOSCF=.F.                  $END
 $MP2    CODE=IMS                                       $END
 $DATA
test GAMESS input file
C1
"""

    file_tail = """ $END
 $I
S   11
1         5.546500E+03           1.560000E-04
2         8.382140E+02           9.860000E-04
3         1.821870E+02           2.792000E-03
4         3.121230E+01          -4.325100E-02
5         1.953140E+01           2.341340E-01
6         8.240990E+00          -7.509430E-01
7         2.194550E+00           8.829680E-01
8         1.109110E+00           4.620610E-01
9         3.746410E-01           2.228600E-02
10        1.770800E-01          -4.353000E-03
11        8.106100E-02           1.102000E-03
S   1
1         3.746410E-01           1.000000E+00
S   1
1         1.770800E-01           1.000000E+00
S   11
1         5.546500E+03          -7.300000E-05
2         8.382140E+02          -5.080000E-04
3         1.821870E+02          -1.158000E-03
4         3.121230E+01           1.219300E-02
5         1.953140E+01          -8.785400E-02
6         8.240990E+00           3.382000E-01
7         2.194550E+00          -5.765500E-01
8         1.109110E+00          -4.092980E-01
9         3.746410E-01           5.674590E-01
10        1.770800E-01           6.124890E-01
11        8.106100E-02           1.432310E-01
S   1
1         1.001000E-01           1.000000E+00
P   9
1         1.889880E+02           5.850000E-04
2         2.128680E+01           3.692300E-02
3         1.003960E+01          -2.353240E-01
4         3.451800E+00           3.414830E-01
5         1.974560E+00           5.347880E-01
6         1.024200E+00           2.651410E-01
7         4.494370E-01           2.578700E-02
8         1.866480E-01           5.220000E-04
9         7.348100E-02           6.060000E-04
P   1
1         5.981000E-01           1.000000E+00
P   9
1         1.889880E+02          -2.560000E-04
2         2.128680E+01          -1.168200E-02
3         1.003960E+01           8.319200E-02
4         3.451800E+00          -1.569700E-01
5         1.974560E+00          -2.245180E-01
6         1.024200E+00          -1.144510E-01
7         4.494370E-01           3.753560E-01
8         1.866480E-01           5.751360E-01
9         7.348100E-02           2.459170E-01
P   1
1         9.618000E-02           1.000000E+00
D   9
1         1.326620E+02           5.720000E-04
2         3.760540E+01           4.402000E-03
3         1.038910E+01          -4.092200E-02
4         6.490170E+00           9.966100E-02
5         3.454510E+00           3.226630E-01
6         1.844130E+00           4.003430E-01
7         9.624780E-01           2.683060E-01
8         4.728530E-01           8.484700E-02
9         1.932000E-01           7.632000E-03
D   1
1         4.728530E-01           1.000000E+00
D   1
1         1.932000E-01           1.000000E+00
F   1
1         4.064000E-01           1.000000E+00

 $END
 $ECP
I-ECP GEN    28    4
1     ----- g-ul potential -----
      0.00000000      2       1.00000000
3     ----- s-g potential -----
     49.98964900      2      40.03337600
    281.00655600      2      17.30057600
     61.41673900      2       8.85172000
4     ----- p-g potential -----
     67.41623900      2      15.72014100
    134.80769600      2      15.20822200
     14.56654800      2       8.29418600
     28.96842200      2       7.75394900
4     ----- d-g potential -----
     35.53875600      2      13.81775100
     53.33975900      2      13.58780500
      9.71646600      2       6.94763000
     14.97750000      2       6.96009900
4     ----- f-g potential -----
    -20.17661800      2      18.52295000
    -26.08807700      2      18.25103500
     -0.22043400      2       7.55790100
     -0.22164600      2       7.59740400
I-ECP
 $END

"""
    assert head == file_head
    assert tail == file_tail


def test_reproduce_from_gau_head_tail():
    gau_file = Path(__file__).parent / "files" / "CH2_opt.com"
    mol = Molecule.from_file(gau_file)
    head, tail = extract_head_tail_from_gau_inp(gau_file)
    new_gau_file = Path(__file__).parent / "files" / "CH2_opt_new.com"
    mol.write_to_gaussian_input(filepath=new_gau_file, head=head, tail=tail)

    gau_file_content = gau_file.read_text()
    new_gau_file_content = new_gau_file.read_text()
    try:
        assert gau_file_content == new_gau_file_content
    finally:
        new_gau_file.unlink()


def test_reproduce_from_gms_head_tail():
    gms_file = Path(__file__).parent / "files" / "I2_mp2.inp"
    mol = Molecule.from_file(gms_file)
    head, tail = extract_head_tail_from_gms_inp(gms_file)
    new_gms_file = Path(__file__).parent / "files" / "I2_mp2_new.inp"
    mol.write_to_gamess_input(filepath=new_gms_file, head=head, tail=tail)

    gms_file_content = gms_file.read_text()
    new_gms_file_content = new_gms_file.read_text()
    try:
        assert gms_file_content == new_gms_file_content
    finally:
        new_gms_file.unlink()
