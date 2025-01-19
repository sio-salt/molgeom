import sys
import re
from pathlib import Path
from collections import deque

from molgeom.atom import Atom
from molgeom.molecule import Molecule
from molgeom.parsers.parser_tools import (
    is_valid_gms_xyz_line,
    remove_trailing_empty_lines,
    validate_filepath,
)


def from_gms_inp_str(content: str) -> Molecule:
    mol = Molecule()

    lines = deque(remove_trailing_empty_lines(content.strip().split("\n")))

    # replace tabs, non-breaking spaces, and multi spaces with single space
    for i in range(len(lines)):
        lines[i] = re.sub(r"[\s\t\xa0]+", " ", lines[i])

    # input description section
    input_descriptions = []
    while lines and not lines[0].strip().upper().startswith("$DATA"):
        if lines[0].strip().startswith("!"):
            lines.popleft()
            continue
        input_descriptions.append(lines.popleft().strip())

    # $DATA group
    lines.popleft()
    if not lines[0].strip() or not lines[1].strip():
        print("inp_parser DATA group : invalid file format \n" + f"{lines[0]}\n" + f"{lines[1]}")
        sys.exit(1)
    _title = lines.popleft().strip()
    _group_naxis = lines.popleft().strip()

    # atom cartesian coordinates
    while lines and lines[0].strip() != "$END":
        if not is_valid_gms_xyz_line(lines[0]):
            print("inp_parser atom cartesian coords : " + f"invalid file format \n{lines[0]}")
            sys.exit(1)
        data = lines.popleft().strip().split()
        atom = Atom(symbol=data[0], x=float(data[2]), y=float(data[3]), z=float(data[4]))
        mol.add_atom(atom)

    return mol


# GAMESS input file parser
def gms_inp_parser(filepath: str | Path) -> Molecule:
    filepath = validate_filepath(filepath)
    with open(filepath, "r") as file:
        content = file.read()
    mol = from_gms_inp_str(content)
    mol.name = filepath.stem
    return mol
