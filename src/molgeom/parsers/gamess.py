import re
from pathlib import Path
from collections import deque

from molgeom.atom import Atom
from molgeom.molecule import Molecule
from molgeom.parsers.parser_tools import (
    is_valid_gms_xyz_line,
    remove_trailing_empty_lines,
    validate_filepath,
    zopen,
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
        raise ValueError(
            "inp_parser DATA group : invalid file format \n"
            + f"{lines[0]}\n"
            + f"{lines[1]}"
        )
    _title = lines.popleft().strip()
    _group_naxis = lines.popleft().strip()

    # atom cartesian coordinates
    while lines and lines[0].strip() != "$END":
        if not is_valid_gms_xyz_line(lines[0]):
            raise ValueError(
                "inp_parser atom cartesian coords :"
                + f"invalid file format \n{lines[0]}"
            )
        data = lines.popleft().strip().split()
        atom = Atom(
            symbol=data[0], x=float(data[2]), y=float(data[3]), z=float(data[4])
        )
        mol.add_atom(atom)

    return mol


def extract_head_tail_from_gms_inp(filepath: str | Path) -> tuple[str, str]:
    file_head = []
    file_tail = []
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"{filepath} do not exist")
    if not filepath.is_file():
        raise ValueError(f"{filepath} is not a file")
    with zopen(filepath, mode="rt", encoding="utf-8") as file:
        while True:
            line = file.readline()
            if line.strip().upper().startswith("$DATA"):
                file_head.append(line)  # $DATA
                line = file.readline()
                file_head.append(line)  # comment line
                line = file.readline()
                if not line.strip().upper().startswith("C1"):
                    raise ValueError(f"Symmetry group {line.strip()} is not supported")
                file_head.append(line)  # symmetry group (e.g. C1)
                break
            file_head.append(line)
        while True:
            line = file.readline()
            if line.strip().upper().startswith("$END"):
                file_tail.append(line)
                break
        while True:
            line = file.readline()
            if not line:
                break
            file_tail.append(line)
    return "".join(file_head), "".join(file_tail)


def gms_inp_parser(filepath: str | Path) -> Molecule:
    filepath = validate_filepath(filepath)
    with zopen(filepath, mode="rt", encoding="utf-8") as file:
        content = file.read()
    mol = from_gms_inp_str(content)
    mol.name = filepath.stem

    return mol
