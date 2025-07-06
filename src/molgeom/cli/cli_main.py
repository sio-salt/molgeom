import sys
from pathlib import Path

import click
from molgeom import Molecule, Vec3, read_file, poscar_parser

if "poscar2xyz" in sys.argv:
    if "--" not in sys.argv:
        idx = sys.argv.index("poscar2xyz") + 2
        if len(sys.argv) > idx:
            sys.argv.insert(idx, "--")


def validate_files(ctx, param, value):
    """Validate existence of input files."""
    files = []
    for file in value:
        path = Path(file.strip()).expanduser().resolve(strict=True)
        if not path.exists():
            raise click.BadParameter(f"File {file} does not exist")
        files.append(str(path))
    return files


@click.group()
def cli():
    """Molecular geometry manipulation tools."""
    pass


@cli.command()
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True),
    required=True,
    callback=validate_files,
)
def center(files):
    """Calculate center of mass for molecular structures."""
    for filepath in files:
        click.echo(f"\n{filepath}")
        mol = read_file(filepath)
        com = mol.center_of_mass()
        click.echo(str(com))


@cli.command()
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True),
    required=True,
    callback=validate_files,
)
def nuclrep(files):
    """Calculate nuclear repulsion energy."""
    for filepath in files:
        click.echo(f"\n{filepath}")
        mol = read_file(filepath)
        energy = mol.nuclear_repulsion()
        click.echo(f" {energy}")


@cli.command()
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True),
    required=True,
    callback=validate_files,
)
@click.option("-t", "--tol", type=float, help="Bond length tolerance (default: 0.15)")
def bonds(files, tol):
    """Analyze molecular bonds with optional tolerance."""
    tolerance = tol
    if tolerance is None:
        tolerance = click.prompt(
            "Enter tolerance of bond length", default=0.15, type=float
        )

    for filepath in files:
        click.echo(f"\n{filepath}")
        mol = read_file(filepath)
        bonds = mol.get_bonds(tolerance)
        click.echo(f"Number of bonds: {len(bonds)}")
        click.echo(f"Bond tolerance: {tolerance}")
        click.echo(
            "label   bond length (Angstrom)       "
            + "atom1                                                                          -     atom2"
        )
        for label, bond_dict in enumerate(bonds):
            i, j = bond_dict["pair"]
            ai, aj = mol[(i, j)]
            dist_angst = ai.distance_to(aj)
            click.echo(
                f"{label:3d}          {dist_angst:.9f}           {i + 1:3d}   {ai}    -   {j + 1:3d}   {aj}"
            )


def translate_molecule(mol):
    """Translate molecule from point A to B."""
    click.echo("\n--- Translation ---")
    pA = click.prompt("Coordinate A", type=str)
    pB = click.prompt("Coordinate B", type=str)
    try:
        pA = Vec3(*map(float, pA.split()))
        pB = Vec3(*map(float, pB.split()))
        trans_vec = pB - pA
        mol.translate(trans_vec)
        return mol
    except (ValueError, TypeError) as e:
        raise click.BadParameter(f"Invalid coordinates: {e}")


def mirror_molecule(mol):
    """Mirror molecule by plane defined by three points."""
    click.echo("\n--- Reflection ---")
    try:
        p1 = Vec3(*map(float, click.prompt("Point 1", type=str).split()))
        p2 = Vec3(*map(float, click.prompt("Point 2", type=str).split()))
        p3 = Vec3(*map(float, click.prompt("Point 3", type=str).split()))
        mol.mirror_by_plane(p1, p2, p3)
        return mol
    except (ValueError, TypeError) as e:
        raise click.BadParameter(f"Invalid points: {e}")


def rotate_molecule(mol):
    """Rotate molecule around axis by angle."""
    click.echo("\n--- Rotation ---")
    try:
        p1 = Vec3(*map(float, click.prompt("Axis point 1", type=str).split()))
        p2 = Vec3(*map(float, click.prompt("Axis point 2", type=str).split()))
        angle = click.prompt("Rotation angle (degrees)", type=float)
        mol.rotate_by_axis(p1, p2, angle)
        return mol
    except (ValueError, TypeError) as e:
        raise click.BadParameter(f"Invalid input: {e}")


OP_FUNCS = {
    "translate": translate_molecule,
    "reflect": mirror_molecule,
    "rotate": rotate_molecule,
}


@cli.command()
@click.argument("file", type=click.Path(exists=True))
@click.option(
    "--operation",
    "-op",
    type=click.Choice(list(OP_FUNCS.keys())),
    help="Geometry modification operation",
)
def modify(file, operation):
    """Modify molecular geometry interactively."""
    mol = read_file(file)
    click.echo(f"\n{file}")

    if operation:
        mol = OP_FUNCS[operation](mol)
    else:
        # Get operation order
        operations = click.prompt(
            "Enter operation order (0: skip, 1-3: order) for [translate reflect rotate]",
            type=str,
        )
        try:
            orders = [int(x) for x in operations.split()]
            if len(orders) != 3 or not all(0 <= x <= 3 for x in orders):
                raise click.BadParameter("Invalid operation order")

            op_order = dict(zip(OP_FUNCS.keys(), orders))
            for i in range(1, 4):
                for op, order in op_order.items():
                    if order == i:
                        mol = OP_FUNCS[op](mol)
        except ValueError:
            raise click.BadParameter("Invalid operation order format")

    click.echo("\nFinal molecule geometry:\n")
    click.echo(mol.to_xyz_str())


@cli.command()
@click.argument("file", type=click.Path(exists=True))
@click.argument("cell_range", nargs=6, type=int)
def poscar2xyz(file, cell_range):
    """
    Convert POSCAR to XYZ format with cell repetition.

    Provide 6 integers: a_min a_max b_min b_max c_min c_max

    Example:
        # 3x3x3 with original cell at center:
        molgeom poscar2xyz POSCAR -1 2 -1 2 -1 2
    """
    mol = poscar_parser(file)
    mol.replicate(
        list(cell_range[0:2]),
        list(cell_range[2:4]),
        list(cell_range[4:6]),
    )
    click.echo(len(mol))
    click.echo(f"{Path(file).stem}, rep: {cell_range}, {str(mol)}")
    click.echo(mol.to_xyz_str())


@cli.command()
@click.argument("file", type=click.Path(exists=True))
def split(file):
    """Split molecular clusters."""
    click.echo(f"\n{file}")
    mol = read_file(file)
    clusters = mol.get_clusters()
    click.echo(f"clusters: {clusters}")
    for cluster in clusters:
        if len(cluster) >= 3:
            click.echo(cluster)
            click.echo(cluster.to_xyz_str())
            click.echo()


@cli.command()
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True),
    required=True,
    callback=validate_files,
)
def view(files):
    """View molecular structure(s) in your browser."""
    mols = []
    for file_path in files:
        try:
            mol = Molecule.from_file(file_path)
            mols.append(mol)
        except Exception as e:
            click.echo(f"Error reading {file_path}: {str(e)}", err=True)
            continue

    if not mols:
        click.echo("No valid molecular geometry files were provided.", err=True)
        return

    Molecule.view_mols(mols)


@cli.command()
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True),
    required=True,
    callback=validate_files,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output file path for combined XYZ",
)
def onexyz(files, output):
    """Combine molecular structures into one XYZ file."""
    mols = []
    for file_path in files:
        try:
            mol = Molecule.from_file(file_path)
            mols.append(mol)
        except Exception as e:
            click.echo(f"Error reading {file_path}: {str(e)}", err=True)
            continue

    if not mols:
        click.echo("No valid molecular geometry files were provided.", err=True)
        return

    Molecule.write_to_xyzs(mols=mols, filepath=output)
    click.echo(f"Combined {len(mols)} molecules into one XYZ file: {output}")
