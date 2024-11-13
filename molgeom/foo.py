from __future__ import annotations
import json
from pathlib import Path

with open(Path(__file__).absolute().parent / "periodic_table.json", encoding="utf-8") as ptable_json:
    ptable = json.load(ptable_json)

C = ptable["C"]

print(C["Atomic mass"])
print(json.dumps(C, indent=4, ensure_ascii=False))
