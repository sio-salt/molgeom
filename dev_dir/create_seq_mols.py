import subprocess
from pathlib import Path
from molgeom import Molecule, Vec3

repo_root = Path(
    subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE)
    .stdout.decode()
    .strip()
)
files_path = repo_root / "tests/files"
output_path = Path(__file__).parent / "seq_mols"
h2o = Molecule.from_file(files_path / "H2O.xyz")  # on YZ plane

rots = []
for i in range(0, 360, 10):
    # rotate around x axis by 10 degrees
    copied = h2o.copy()
    copied.rotate_by_axis(axis_point1=[0, 0, 0], axis_point2=[1, 0, 0], deg=i)
    rots.append(copied)

for i, mol in enumerate(rots):
    mol.write_to_xyz(output_path / f"rot_h2o_{i:03d}.xyz")
