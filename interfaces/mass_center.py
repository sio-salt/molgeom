import sys

from molgeom import read_file


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_path1> <file_path2> ...")
        sys.exit(1)

    print()
    filepaths = sys.argv[1:]
    for filepath in filepaths:
        print(filepath)
        mole = read_file(filepath)
        com = mole.center_of_mass()
        print(com)
        print()


if __name__ == "__main__":
    main()
