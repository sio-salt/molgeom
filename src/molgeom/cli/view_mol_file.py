import argparse
from pathlib import Path
from molgeom import Molecule


def main():
    parser = argparse.ArgumentParser(description="View molecular geometry file(s) using 3Dmol.js")
    parser.add_argument(
        "files",
        nargs="+",
        type=str,
        help="Path to molecular geometry file(s) (xyz, com, gjf, cif, POSCAR etc.)",
    )

    args = parser.parse_args()
    molecules = []

    for file_path in args.files:
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: File {file_path} does not exist, skipping...")
            continue
        try:
            mol = Molecule.from_file(file_path)
            molecules.append(mol)
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            continue

    if not molecules:
        print("No valid molecular geometry files were provided.")
        return

    Molecule.view_mols(molecules)


if __name__ == "__main__":
    main()
