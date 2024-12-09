import sys
from molgeom import parse_file


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <xyz_file_path> ...")
        sys.exit(1)

    print()
    filepath = sys.argv[1]
    print(filepath)
    mole = parse_file(filepath)
    clusters = mole.get_clusters()
    print(f"{clusters=}")
    for cluster in clusters:
        if len(cluster) >= 3:
            print(cluster)
            print(cluster.to_xyz())


if __name__ == "__main__":
    main()
