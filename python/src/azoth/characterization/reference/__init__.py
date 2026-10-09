"""The pure-Python reference implementation for the characterisation namespace.

Same contract as the other namespaces' references: always importable, with or without the
compiled extension, and written to mirror the Rust line for line so a reviewer can read the two
side by side against the correlations the spec cites.
"""

from __future__ import annotations

# Imported first so the submodule is bound before the kernels below reach for it, and as a
# redundant alias because that is what `ruff` needs to see the import as deliberate.
from azoth.characterization.reference import _unported as _unported
from azoth.characterization.reference.lumping import lumping
from azoth.characterization.reference.pedersen_plus_split import pedersen_plus_split
from azoth.characterization.reference.tbp_closure import tbp_closure
from azoth.characterization.reference.tbp_cut_properties import tbp_cut_properties
from azoth.characterization.reference.tbp_density import tbp_density
from azoth.characterization.reference.whitson_gamma_split import whitson_gamma_split

__all__ = [
    "lumping",
    "pedersen_plus_split",
    "tbp_closure",
    "tbp_cut_properties",
    "tbp_density",
    "whitson_gamma_split",
]
