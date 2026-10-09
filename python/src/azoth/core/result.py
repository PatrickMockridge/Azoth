"""The result shapes: the generated dataclasses, and the hand-written parts they are built from.

The 194 result classes are in :mod:`azoth.core.result_gen`, emitted by
``tools/gen_python_result.py`` from three declarations - ``CalcResult::FIELDS`` for the fields and
their order, the transport's Rust type beside the spec's declared unit for each annotation, and the
spec's own ``description`` for each doc. The enums those classes annotate, the ``_HasWarnings``
mixin they inherit, and the fitting record that is not a result are in
:mod:`azoth.core.result_base`, where they are written by hand.

**The split is what breaks a cycle.** A generated module cannot import the file that re-exports it,
so the base has to live somewhere the generated file can reach without going through this one - the
arrangement :mod:`azoth.core.bridge_support` makes for the Rust bridge.

This module stays the import site for both, so ``from azoth.core.result import FlowRegime`` and
``from azoth.core.result import HaalandResult`` are unchanged.
"""

from __future__ import annotations

from azoth.core.result_base import LAMINAR_MAX as LAMINAR_MAX
from azoth.core.result_base import TURBULENT_MIN as TURBULENT_MIN
from azoth.core.result_base import FlowRegime as FlowRegime
from azoth.core.result_base import HenryStatus as HenryStatus
from azoth.core.result_base import HydrateStructure as HydrateStructure
from azoth.core.result_base import KComponent as KComponent
from azoth.core.result_base import Phase as Phase
from azoth.core.result_base import PlusModel as PlusModel
from azoth.core.result_base import RootStructure as RootStructure
from azoth.core.result_base import StabilityVerdict as StabilityVerdict
from azoth.core.result_base import TpMultiflashSeed as TpMultiflashSeed
from azoth.core.result_base import _HasWarnings as _HasWarnings
from azoth.core.result_gen import *  # noqa: F403
