import importlib.util
import os
import sys

from molgeom import read_file

if importlib.util.find_spec("readline"):
    import readline  # noqa


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)
    else:
        filepath = sys.argv[1]
        if not os.path.isfile(filepath):
            print(f"File not found: {filepath}")
            sys.exit(1)

    print()
    while True:
        try:
            print(
                'Enter tolerance of bond length ( e.g. 0.15  or "default": [default tol = 0.15 (Å)] )'
            )
            inp_str = input(" >>> ").strip()
            if inp_str.lower() == "default":
                tolerance = 0.15
            else:
                tolerance = float(inp_str)

            break
        except ValueError:
            print("Please Enter valid float value")
        # lower_bound = float(input("Enter the lower bond threshold: "))
        # upper_bound = float(input("Enter the upper bond threshold: "))

    print()
    filepath = sys.argv[1]
    print(filepath)
    output = []
    mole = read_file(filepath)
    bonds = mole.get_bonds(tolerance)
    output.append(f"Number of bonds: {len(bonds)}")
    output.append(f"bond tolerance: {tolerance}")
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

    # output_file = "bonds_output.txt"
    # with open(output_file, "w") as f:
    #     f.write("\n".join(output))
    #
    # print(f"Output written to {output_file}")
    # print()


if __name__ == "__main__":
    main()
