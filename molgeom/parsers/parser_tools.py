import re

from molgeom.data.consts import ATOMIC_MASSES, SPECIAL_ELEMENTS

# atom_xyz_regex = re.compile(r"(\w\w?)(\s+[-+]?\d*\.\d+){3}")
symbol_xyz_regex = re.compile(r"(\w\w?)(\s+[-+]?(?:\d+|\d*\.\d+)){3}")
symbol_mass_xyz_regex = re.compile(r"(\w\w?)(\s+\d+\.\d+)(\s+[-+]?\d*\.\d+){3}")


def remove_trailing_empty_lines(lines: list[str]):
    while lines and not lines[-1].strip():
        lines.pop()
    return lines


def is_valid_xyz_line(line: str) -> bool:
    data = line.strip().split()
    if not data:
        return False
    if len(data) != 4:
        return False
    if not re.fullmatch(symbol_xyz_regex, line.strip()):
        return False
    if not (data[0] in ATOMIC_MASSES or data[0] in SPECIAL_ELEMENTS):
        return False
    return True


def is_valid_gms_xyz_line(line: str) -> bool:
    data = line.strip().split()
    if not data:
        return False
    if len(data) != 5:
        return False
    if not re.fullmatch(symbol_mass_xyz_regex, line.strip()):
        return False
    if data[0] not in ATOMIC_MASSES:
        return False
    return True
