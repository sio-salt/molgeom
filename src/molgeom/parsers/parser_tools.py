import re
import gzip
import bz2
import lzma
import warnings
from pathlib import Path


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


def validate_filepath(filepath: str | Path) -> Path:
    filepath = Path(str(filepath).strip()).expanduser().resolve(strict=True)
    if not filepath.exists():
        raise FileNotFoundError(f"{filepath} do not exist")
    if not filepath.is_file():
        raise ValueError(f"{filepath} is not a file")
    return filepath


def zopen(file: str | Path, mode=None, **kwargs):
    """
    Open a file with transparent support for compressed files (.gz, .bz2, .xz/.lzma).
    This version uses pathlib for path handling.

    Args:
        filename (str or Path): The file path.
        mode (str, optional): Mode to open the file.
            Must explicitly include 't' (text) or 'b' (binary).
            If not specified, defaults to 'r' with a warning.
        **kwargs: Additional keyword arguments for the underlying open functions.

    Returns:
        A file-like object.
    """
    filepath = Path(file)

    if mode is None:
        warnings.warn(
            "Explicit mode not specified. Defaulting to 'r'. "
            "Please specify mode including 't' (text) or 'b' (binary).",
            FutureWarning,
            stacklevel=2,
        )
        mode = "r"

    if "t" not in mode and "b" not in mode:
        warnings.warn(
            "Mode should explicitly include 't' (text) or 'b' (binary).",
            FutureWarning,
            stacklevel=2,
        )

    ext = filepath.suffix.lower()
    if ext == ".gz":
        return gzip.open(str(filepath), mode, **kwargs)
    elif ext == ".bz2":
        return bz2.open(str(filepath), mode, **kwargs)
    elif ext in (".xz", ".lzma"):
        return lzma.open(str(filepath), mode, **kwargs)
    else:
        return filepath.open(mode, **kwargs)
