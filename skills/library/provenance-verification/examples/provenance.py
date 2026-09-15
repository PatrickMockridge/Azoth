"""Where a result's provenance lives: NOTICE and the warnings on a result."""

from pathlib import Path

import azoth

r = azoth.hydraulics.darcy_weisbach(
    0.02,
    azoth.ureg.Quantity(100.0, "m"),
    azoth.ureg.Quantity(0.1, "m"),
    azoth.ureg.Quantity(998.0, "kg/m**3"),
    azoth.ureg.Quantity(1.5, "m/s"),
)
print("is_clean:", r.is_clean)
print("warnings:", [str(w.code) for w in r.warnings])

notice = Path(__file__).resolve().parents[4] / "NOTICE"
print("NOTICE lines:", len(notice.read_text(encoding="utf-8").splitlines()))
