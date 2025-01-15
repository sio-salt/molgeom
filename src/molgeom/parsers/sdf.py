from pathlib import Path

from molgeom.molecule import Molecule
from molgeom.parsers.mol import from_mol_str


def from_sdf_str(content: str) -> list[Molecule]:
    """Parse SDF format string content and return list of Molecule objects.

    SDF format contains multiple MOL format molecules separated by $$$$
    """
    molecules = []
    # Split content at $$$$ to get individual MOL blocks
    mol_blocks = content.split("$$$$")

    for block in mol_blocks:
        block = block.strip()
        if block:  # Skip empty blocks
            try:
                mol = from_mol_str(block)
                molecules.append(mol)
            except ValueError as e:
                raise ValueError(f"Error parsing MOL block: {e}")

    if not molecules:
        raise ValueError("No valid molecules found in SDF content")

    return molecules


def sdf_parser(filepath: str | Path) -> list[Molecule]:
    """Read SDF file and return list of Molecule objects."""
    filepath = str(filepath)
    with open(filepath, "r") as f:
        content = f.read()
    return from_sdf_str(content)
