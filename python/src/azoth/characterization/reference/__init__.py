"""The pure-Python reference implementation for the characterisation namespace.

Same contract as the other namespaces' references: always importable, with or without the
compiled extension, and written to mirror the Rust line for line so a reviewer can read the two
side by side against the correlations the spec cites.
"""

from __future__ import annotations

from azoth.characterization.reference.tbp_cut_properties import tbp_cut_properties

__all__ = [
    "tbp_cut_properties",
]
