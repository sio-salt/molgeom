import click
from pathlib import Path
import sys

from molgeom import read_file, poscar_parser, Molecule


def validate_files(ctx, param, value):
    """ファイルの存在確認を行うコールバック関数"""
    files = []
    for file in value:
        path = Path(file)
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
    "files", nargs=-1, type=click.Path(exists=True), required=True, callback=validate_files
)
def center(files):
    """Calculate center of mass for molecular structures."""
    for filepath in files:
        click.echo(f"\n{filepath}")
        mole = read_file(filepath)
        com = mole.center_of_mass()
        click.echo(str(com))


@cli.command()
@click.argument(
    "files", nargs=-1, type=click.Path(exists=True), required=True, callback=validate_files
)
def nuclrep(files):
    """Calculate nuclear repulsion energy."""
    for filepath in files:
        click.echo(f"\n{filepath}")
        mole = read_file(filepath)
        energy = mole.nuclear_repulsion()
        click.echo(f" {energy}")


@cli.command()
@click.argument(
    "files", nargs=-1, type=click.Path(exists=True), required=True, callback=validate_files
)
@click.option("-t", "--tol", type=float, help="Bond length tolerance (default: 0.15)")
def bonds(files, tol):
    """Analyze molecular bonds with optional tolerance."""
    tolerance = tol
    if tolerance is None:
        tolerance = click.prompt("Enter tolerance of bond length", default=0.15, type=float)

    for filepath in files:
        click.echo(f"\n{filepath}")
        mole = read_file(filepath)
        bonds = mole.get_bonds(tolerance)
        click.echo(f"Number of bonds: {len(bonds)}")
        click.echo(f"Bond tolerance: {tolerance}")
        click.echo(
            "label   bond length (Angstrom)       "
            + "atom1                                                                          -     atom2"
        )
        for label, bond_dict in enumerate(bonds):
            i, j = bond_dict["pair"]
            ai, aj = mole[(i, j)]
            dist_angst = ai.distance_to(aj)
            click.echo(
                f"{label:3d}          {dist_angst:.9f}           {i+1:3d}   {ai}    -   {j+1:3d}   {aj}"
            )


@cli.command()
@click.argument("file", type=click.Path(exists=True))
@click.option(
    "--operation",
    "-op",
    type=click.Choice(["translate", "reflect", "rotate"]),
    help="Geometry modification operation",
)
def modify(file, operation):
    """Modify molecular geometry interactively."""
    from cli.geom_modifier import (
        translate_molecule,
        mirror_molecule,
        rotate_molecule,
        get_operation_order,
    )

    click.echo(f"\n{file}\n")
    mole = read_file(file)

    if operation:
        op_funcs = {
            "translate": translate_molecule,
            "reflect": mirror_molecule,
            "rotate": rotate_molecule,
        }
        mole = op_funcs[operation](mole)
    else:
        operation_order = get_operation_order()
        if all(v == 0 for v in operation_order.values()):
            sys.exit(0)

        for i in range(1, 4):
            for op, order in operation_order.items():
                if order == i:
                    op_funcs = {
                        "translate": translate_molecule,
                        "reflect": mirror_molecule,
                        "rotate": rotate_molecule,
                    }
                    mole = op_funcs[op](mole)

    click.echo("\n final molecule geometry \n")
    click.echo(mole.to_xyz())


@cli.command()
@click.argument("file", type=click.Path(exists=True))
@click.argument("cell_range", nargs=6, type=int)
def poscar2xyz(file, cell_range):
    """Convert POSCAR to XYZ format with cell repetition.

    Cell range should be specified as: a_min a_max b_min b_max c_min c_max
    """
    mole = poscar_parser(file)
    mole.replicate(
        cell_range[0:2],
        cell_range[2:4],
        cell_range[4:6],
    )
    click.echo(mole.to_xyz())


@cli.command()
@click.argument("file", type=click.Path(exists=True))
def split(file):
    """Split molecular clusters."""
    click.echo(f"\n{file}")
    mole = read_file(file)
    clusters = mole.get_clusters()
    click.echo(f"clusters: {clusters}")
    for cluster in clusters:
        if len(cluster) >= 3:
            click.echo(cluster)
            click.echo(cluster.to_xyz())
            click.echo()


@cli.command()
@click.argument(
    "files", nargs=-1, type=click.Path(exists=True), required=True, callback=validate_files
)
def view(files):
    """View molecular structure(s)."""
    molecules = []
    for file_path in files:
        try:
            mol = Molecule.from_file(file_path)
            molecules.append(mol)
        except Exception as e:
            click.echo(f"Error reading {file_path}: {str(e)}", err=True)
            continue

    if not molecules:
        click.echo("No valid molecular geometry files were provided.", err=True)
        return

    Molecule.view_mols(molecules)


if __name__ == "__main__":
    cli()
