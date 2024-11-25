import sys
from molgeom import parse_file


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <xyz_file_path> ...")
        sys.exit(1)

    print()
    filepath = sys.argv[1]
    print(filepath)
    mole = parse_file(filepath, "r")
    clusters = mole.get_bond_clusters()
    print(f"{clusters=}")
    for cluster in clusters:
        print(cluster.to_xyz())
        print()


if __name__ == "__main__":
    main()