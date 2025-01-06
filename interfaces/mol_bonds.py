import importlib.util
import os
import sys
import argparse
from molgeom import read_file

if importlib.util.find_spec("readline"):
    import readline  # noqa


def main():
    parser = argparse.ArgumentParser(
        description="Calculate molecular bonds with optional tolerance."
    )
    parser.add_argument(
        "files",
        metavar="FILE",
        nargs="+",
        help="Path(s) to the molecular structure file(s).",
    )
    parser.add_argument(
        "-t",
        "--tol",
        type=float,
        default=None,
        help="Tolerance of bond length (e.g., 0.15). If not provided, you will be prompted to enter it.",
    )
    args = parser.parse_args()

    for filepath in args.files:
        if not os.path.isfile(filepath):
            print(f"File not found: {filepath}")
            sys.exit(1)

    if args.tol is not None:
        tolerance = args.tol
    else:
        print()
        while True:
            try:
                print(
                    'Enter tolerance of bond length (e.g., 0.15 or "default": [default tol = 0.15 (Å)] )'
                )
                inp_str = input(" >>> ").strip()
                if inp_str.lower() == "default":
                    tolerance = 0.15
                else:
                    tolerance = float(inp_str)
                break
            except ValueError:
                print("Please enter a valid float value.")

    print()
    for filepath in args.files:
        print(filepath)
        output = []
        mole = read_file(filepath)
        bonds = mole.get_bonds(tolerance)
        output.append(f"Number of bonds: {len(bonds)}")
        output.append(f"Bond tolerance: {tolerance}")
        output.append(
            "label   bond length (Angstrom)       "
            + "atom1                                                                          -     atom2"
        )
        for label, bond_dict in enumerate(bonds):
            i, j = bond_dict["pair"]
            ai, aj = mole[(i, j)]
            dist_angst = ai.distance_to(aj)
            output.append(
                f"{label:3d}          {dist_angst:.9f}           {i+1:3d}   {ai}    -   {j+1:3d}   {aj}"
            )
        print("\n".join(output))
        print()

    # Uncomment to write output to a file
    # output_file = "bonds_output.txt"
    # with open(output_file, "w") as f:
    #     f.write("\n".join(output))
    # print(f"Output written to {output_file}")


if __name__ == "__main__":
    main()
