"""Components and mixtures: the inputs a flash takes.

Two small data holders rather than arithmetic. The arithmetic lives in
``azoth.eos.reference.pt_flash`` (and its Rust counterpart), and the shape of a
mixture's input is the same for both backends - so it is defined once, here, where
neither implementation owns it.

# No component names, deliberately

:class:`Component` carries ``Tc``, ``Pc`` and ``omega`` and **no name**: a name would
let a *calculation* resolve it, and a flash that looks up its own inputs is one whose
answer depends on a file the caller never mentioned. So the lookup is a step the
caller takes - ``azoth.eos.components.component("methane")`` returns one of these, and
everything downstream still takes the numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from azoth.core.errors import InvalidInputError
from azoth.eos.cubic import PR, Cubic

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from azoth.core.units import Q


@dataclass(frozen=True, slots=True)
class Component:
    """One component's critical constants.

    All three are the caller's. A table of them is a licensing and provenance
    problem as much as a data problem, and shipping one is a decision worth making
    on its own rather than in passing - so a :class:`Component` is the only place
    they exist.

    Attributes:
        Tc: critical temperature.
        Pc: critical pressure.
        omega: acentric factor.
    """

    Tc: Q
    Pc: Q
    omega: float

    def __post_init__(self) -> None:
        for name in ("Tc", "Pc"):
            value = getattr(self, name).to_base_units().magnitude
            if not value > 0.0:
                raise InvalidInputError(
                    name,
                    f"a critical {name} is a divisor in the reduced variables and is "
                    f"positive by definition, but got {value}",
                )


@dataclass(frozen=True, slots=True)
class Mixture:
    """A set of components and the binary interaction parameters between them.

    ``kij`` defaults to zero everywhere - the ideal-mixture case, and the one a
    caller who has no fitted parameters should be using rather than a guessed
    number. It may be given either as a full ``N x N`` matrix or as a mapping of
    ``(i, j)`` pairs to values, which is the form a literature value arrives in.

    Attributes:
        components: one :class:`Component` per component, in the order every other
            argument to :func:`azoth.eos.pt_flash` is indexed by.
        kij: the interaction parameters. Zero diagonal, and symmetric; see below.
    """

    components: tuple[Component, ...]
    kij: tuple[tuple[float, ...], ...] = field(default=())
    cubic: Cubic = field(default=PR)

    def __post_init__(self) -> None:
        if not self.components:
            raise InvalidInputError("components", "a mixture needs at least one component")
        n = len(self.components)
        if self.kij:
            if len(self.kij) != n or any(len(row) != n for row in self.kij):
                raise InvalidInputError(
                    "kij",
                    f"a mixture of {n} components needs an N x N matrix, but the "
                    f"matrix given is {len(self.kij)} x "
                    f"{len(self.kij[0]) if self.kij else 0}",
                )
            for i in range(n):
                if self.kij[i][i] != 0.0:
                    raise InvalidInputError(
                        "kij",
                        f"kij[{i}][{i}] is {self.kij[i][i]} but the diagonal must be "
                        f"zero: a component does not interact with itself, and a "
                        f"non-zero diagonal silently rescales that component's "
                        f"attraction",
                    )
                for j in range(i + 1, n):
                    if self.kij[i][j] != self.kij[j][i]:
                        raise InvalidInputError(
                            "kij",
                            f"kij[{i}][{j}] is {self.kij[i][j]} but kij[{j}][{i}] is "
                            f"{self.kij[j][i]}; the matrix must be symmetric",
                        )
        else:
            object.__setattr__(self, "kij", tuple(tuple(0.0 for _ in range(n)) for _ in range(n)))

    def __len__(self) -> int:
        """How many components there are."""
        return len(self.components)

    def flattened_kij(self) -> list[float]:
        """The interaction matrix as one row-major list, for the Rust boundary."""
        return [value for row in self.kij for value in row]


def mixture(
    components: Iterable[Component],
    kij: Mapping[tuple[int, int], float] | None = None,
    cubic: Cubic = PR,
) -> Mixture:
    """A :class:`Mixture` from a component list and sparse interaction pairs.

    The spelling a caller with a couple of published ``kij`` values actually wants::

        mixture([methane, butane], kij={(0, 1): 0.05})

    Args:
        components: one :class:`Component` per component, in the order the flash's
            ``z`` will be indexed by.
        kij: interaction parameters keyed by ``(i, j)`` with ``i < j``. Each pair may
            be given in either order. Omitted pairs are zero.
        cubic: the cubic the mixture is evaluated under, from
            :mod:`azoth.eos.cubic`. Defaults to Peng-Robinson.

    Returns:
        The mixture, with a full symmetric matrix built from the pairs.

    Raises:
        InvalidInputError: if the pairs name an index outside the component list.
    """
    resolved = tuple(components)
    n = len(resolved)
    matrix = [[0.0] * n for _ in range(n)]
    for (i, j), value in (kij or {}).items():
        if not (0 <= i < n and 0 <= j < n):
            raise InvalidInputError(
                "kij",
                f"the pair ({i}, {j}) names a component outside the {n} supplied",
            )
        if i == j:
            raise InvalidInputError(
                "kij",
                f"the diagonal entry ({i}, {i}) must be zero; a component does not "
                f"interact with itself",
            )
        matrix[i][j] = value
        matrix[j][i] = value
    return Mixture(
        components=resolved,
        kij=tuple(tuple(row) for row in matrix),
        cubic=cubic,
    )
