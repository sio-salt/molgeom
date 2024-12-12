import sys

from molgeom import parse_file


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <xyz_file_path1> <xyz_file_path2> ...")
        sys.exit(1)

    print()
    filepaths = sys.argv[1:]
    for filepath in filepaths:
        print(filepath)
        mole = parse_file(filepath)
        nuclrep = mole.nuclear_repulsion()
        print(f" {nuclrep}")
        print()


if __name__ == "__main__":
    main()
