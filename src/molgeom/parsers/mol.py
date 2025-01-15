from pathlib import Path
from collections import deque
from molgeom.atom import Atom
from molgeom.molecule import Molecule


def from_mol_str(content: str) -> Molecule:
    """Parse MOL format string content and return Molecule object.

    Handles V2000 format MOL files with atom coordinates and charges.
    Bonds are ignored as per requirements.
    """
    mol = Molecule()
    lines = deque(content.strip().split("\n"))

    # Skip header (3 lines: ID, program info, comment)
    for _ in range(3):
        lines.popleft()

    # Parse counts line
    counts = lines.popleft()
    try:
        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])
    except (ValueError, IndexError):
        raise ValueError(f"Invalid counts line: {counts}")

    # Parse atom block
    for _ in range(n_atoms):
        line = lines.popleft().strip()
        try:
            # V2000 format: coordinates are 10 char fields, followed by atom symbol
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            symbol = line[31:34].strip()

            # Create atom and add it to molecule
            atom = Atom(symbol=symbol, x=x, y=y, z=z)
            mol.add_atom(atom)

        except (ValueError, IndexError) as e:
            raise ValueError(f"Invalid atom line: {line}\nError: {e}")

    return mol


def mol_parser(filepath: str | Path) -> Molecule:
    """Read MOL file and return Molecule object."""
    filepath = str(filepath)
    with open(filepath, "r") as f:
        content = f.read()
    return from_mol_str(content)
