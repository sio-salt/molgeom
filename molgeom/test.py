# from molgeom import parse_file
import json
from pathlib import Path

with open(
    Path(__file__).absolute().parent / "periodic_table.json", encoding="utf-8"
) as ptable_json:
    _pt_data = json.load(ptable_json)

print(_pt_data["C"])

# file_path = "/home/kato/10.git_repos/molgeom/exapmle_inputs/1.input_example.xyz"
# mole = parse_file(file_path, "r")

# print(mole[3])
# print(mole[3:5])
