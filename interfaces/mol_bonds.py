import os
import sys
from molgeom import parse_file


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <xyz_file_path1> <xyz_file_path2> ...")
        sys.exit(1)
    else:
        for filepath in sys.argv[1:]:
            if not os.path.isfile(filepath):
                print(f"File not found: {filepath}")
                sys.exit(1)

    print()
    lower_bound = float(input("Enter the lower bond threshold: "))
    upper_bound = float(input("Enter the upper bond threshold: "))

    print()
    filepaths = sys.argv[1:]
    for filepath in filepaths:
        print(filepath)
        output = []
        mole = parse_file(filepath, "r")
        bonds = mole.bonds(lower_bound, upper_bound)
        output.append(f"Number of bonds: {len(bonds)}")
        output.append(f"bond threshold,   lower: {lower_bound}, upper: {upper_bound}")
        output.append(
            "label   bond length (Angstrom)       atom1                                                                 -     atom2"
        )
        for label, (i, j, dist_angst, atom1, atom2) in enumerate(bonds, start=1):
            output.append(
                f"{label:3d}          {dist_angst:.9f}           {i+1:3d}   {atom1}    -   {j+1:3d}   {atom2}"
            )
        print("\n".join(output))
        print()

    output_file = "bonds_output.txt"
    with open(output_file, "w") as f:
        f.write("\n".join(output))

    print(f"Output written to {output_file}")
    print()


if __name__ == "__main__":
    main()

