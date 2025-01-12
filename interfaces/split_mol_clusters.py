import sys

from molgeom import read_file


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    print()
    filepath = sys.argv[1]
    print(filepath)
    mole = read_file(filepath)
    clusters = mole.get_clusters()
    print(f"{clusters=}")
    for cluster in clusters:
        if len(cluster) >= 3:
            print(cluster)
            print(cluster.to_xyz())
            print()


if __name__ == "__main__":
    main()
