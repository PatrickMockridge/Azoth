"""Show the port pipeline: the generators that turn a spec into generated files."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]

generators = sorted(p.name for p in (ROOT / "tools").glob("gen_*.py"))
specs = sorted((ROOT / "specs" / "calcs").rglob("*.toml"))

print("generators:", ", ".join(generators))
print(f"{len(specs)} calc spec(s) under specs/calcs/")
