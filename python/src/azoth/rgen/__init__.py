"""Emit the documents and on-chain payloads that file an azoth run with r-gen.

A bundle is evidence, not a client: it computes no physics of its own beyond running the
registered calculation on both implementations, and it has no dependency on r-gen. The
r-gen harness reads the bundle, checks it and signs the deploys.

Typical use::

    from azoth.rgen import build

    bundle = build(
        "eos.antoine_vapor_pressure",
        {"A": 5.23243, "B": 891.0098, "C": 332.0975, "D": 0.0, "E": 0.0,
         "form": "pow10", "Tc": 190.56, "Pc": 4_599_000.0, "T": 150.0},
        project="azoth-demo",
        jurisdiction="US-TX",
    )
    bundle.write(Path("artifacts/antoine"))
"""

from __future__ import annotations

from .bundle import Bundle, BundleRefused, build
from .canonical import PayloadError
from .inputs import InputError, UnknownCalc, UnsupportedCalc

__all__ = [
    "Bundle",
    "BundleRefused",
    "InputError",
    "PayloadError",
    "UnknownCalc",
    "UnsupportedCalc",
    "build",
]
