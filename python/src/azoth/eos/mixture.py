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

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from azoth.core.errors import InvalidInputError
from azoth.eos.cubic import PR, Cubic

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

    from azoth.core.units import Q
    from azoth.eos.components import (
        AssociationParameters,
        FurstElectrolyte,
        HuronVidalParameters,
        UnifacUmrpruParameters,
    )
    from azoth.eos.reference._mixture_state import ReducedParameters as ReducedParametersLike


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
        alpha_params: fitted parameters for the alpha correlations that need them,
            in the order the correlation reads them. Empty for a component built
            without a fitted set, which is every caller-supplied component and every
            correlation that needs none.
        association: this substance's association parameters, or ``None`` for one
            that carries no site scheme. Carried whether or not the mixture runs
            them - see :attr:`Mixture.associating`, which is the *model's* decision.
    """

    Tc: Q
    Pc: Q
    omega: float
    molar_mass: Q | None = None
    alpha_params: tuple[float, ...] = ()
    association: AssociationParameters | None = None

    def __post_init__(self) -> None:
        for name in ("Tc", "Pc"):
            value = getattr(self, name).to_base_units().magnitude
            if not value > 0.0:
                raise InvalidInputError(
                    name,
                    f"a critical {name} is a divisor in the reduced variables and is "
                    f"positive by definition, but got {value}",
                )


#: The roles a component can take in the Soreide-Whitson aqueous correlation, as the
#: strings :func:`soreide_whitson_role` returns. `hydrocarbon` is NeqSim's `else` branch, so
#: it is the role of anything the other three names do not claim.
SOREIDE_WHITSON_ROLES: tuple[str, ...] = ("water", "nitrogen", "carbon_dioxide", "hydrocarbon")


def soreide_whitson_role(name: str) -> str:
    """Which role a component's **name** gives it in the Soreide-Whitson correlation.

    NeqSim's own test, on `EosMixingRuleHandler.getkijWhitsonSoreideAqueous`: every role is
    decided by the name, because the correlation was fitted per named gas. `H2O` is water
    beside `water`, `N2` beside `nitrogen`, and `CO2` has no synonym there, so a mixture
    naming its carbon dioxide anything else takes the hydrocarbon branch.
    """
    lowered = name.strip().lower()
    if lowered in ("water", "h2o"):
        return "water"
    if lowered in ("n2", "nitrogen"):
        return "nitrogen"
    if lowered == "co2":
        return "carbon_dioxide"
    return "hydrocarbon"


def soreide_whitson_roles(names: Iterable[str]) -> tuple[str, ...]:
    """One role per component, from the names, in the order given."""
    return tuple(soreide_whitson_role(name) for name in names)


@dataclass(frozen=True, slots=True)
class SoreideWhitsonParameters:
    """The two things the Soreide-Whitson rule carries beside the base matrix.

    The base matrix is the mixture's own :attr:`Mixture.kij`, which for this rule is
    `INTER.csv`'s `KIJWhitsonSoriede` column.

    Attributes:
        roles: one role per component, from :func:`soreide_whitson_roles`. Carried rather
            than re-derived from a name because :class:`Component` carries no name - the
            same reason Rust's rule carries a `Vec<SoreideWhitsonRole>`.
        salinity: the equivalent-NaCl molality of the aqueous phase, in mol/kg water.
    """

    roles: tuple[str, ...]
    salinity: float

    def __post_init__(self) -> None:
        for index, role in enumerate(self.roles):
            if role not in SOREIDE_WHITSON_ROLES:
                raise InvalidInputError(
                    "soreide_whitson",
                    f"role {index} is {role!r}; the roles are {list(SOREIDE_WHITSON_ROLES)}",
                )
        if self.salinity < 0.0:
            raise InvalidInputError(
                "soreide_whitson",
                f"the salinity is {self.salinity} mol/kg, and a molality cannot be negative",
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
        associating: whether a phase model runs the Wertheim association
            contribution over these components. **False by default, and opt-in for
            the reason the components carry their parameters unconditionally**: the
            same methanol and water are an associating fluid under `SystemSrkCPA`
            and a classical one under `SystemNRTL`, so whether a mixture associates
            is the model's decision and not a property of the substances.
    """

    components: tuple[Component, ...]
    kij: tuple[tuple[float, ...], ...] = field(default=())
    cubic: Cubic = field(default=PR)
    alpha: str = field(default="pr")
    associating: bool = field(default=False)
    #: How the components' attraction and covolume combine. ``"classic"`` is the van der
    #: Waals one-fluid rule, which reads :attr:`kij`; ``"umr"`` is the universal mixing
    #: rule, which reads the components' UNIFAC group decomposition instead and does not
    #: read :attr:`kij` at all - so a mixture naming it carries a matrix that nothing
    #: looks at.
    mixing_rule: str = field(default="classic")
    #: The UNIFAC-UMR-PRU tables the ``"umr"`` rule reads, or ``None`` for every mixture
    #: that uses an interaction matrix. Carried on the mixture rather than looked up per
    #: state for the reason Rust's ``MixingRule::Umr`` carries them: the tables are a
    #: property of the fluid, and resolving them inside the solve would put a databank
    #: read on a flash's inner loop.
    umr: UnifacUmrpruParameters | None = field(default=None)
    #: The Soreide-Whitson rule's roles and salinity, or ``None`` for every mixture that
    #: does not name that rule. Its *interaction matrix* is :attr:`kij` - the column the
    #: rule reads as its base - and the correlation replaces the water-gas entries of a
    #: water-rich phase with its own. See :meth:`phase_kij`.
    soreide_whitson: SoreideWhitsonParameters | None = field(default=None)
    #: The Huron-Vidal rule's fitted NRTL matrices, or ``None`` for every mixture that does
    #: not name that rule. A parameterised rule like the two beside it, and for the same
    #: reason: its parameters are a property of the fluid.
    huron_vidal: HuronVidalParameters | None = field(default=None)
    #: The Fürst electrolyte term, or ``None`` for every mixture that is not one. Carried
    #: whole rather than as a flag because its short-range table is built from the *whole*
    #: composition's names.
    furst: FurstElectrolyte | None = field(default=None)

    def __post_init__(self) -> None:
        if not self.components:
            raise InvalidInputError("components", "a mixture needs at least one component")
        n = len(self.components)
        if self.mixing_rule == "huron_vidal" and self.huron_vidal is None:
            raise InvalidInputError(
                "huron_vidal",
                "the Huron-Vidal rule mixes the attraction with a co-volume-weighted NRTL "
                "whose parameters are the fluid's, and none is set",
            )
        if self.mixing_rule == "soreide_whitson":
            if self.soreide_whitson is None:
                raise InvalidInputError(
                    "soreide_whitson",
                    "the Soreide-Whitson rule resolves its water-gas entries from each "
                    "component's role and the brine's salinity, and neither is set",
                )
            if len(self.soreide_whitson.roles) != n:
                raise InvalidInputError(
                    "soreide_whitson",
                    f"the mixture has {n} components and the rule carries "
                    f"{len(self.soreide_whitson.roles)} roles",
                )
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

    def phase_kij(
        self, reduced: ReducedParametersLike, x: Sequence[float]
    ) -> tuple[tuple[float, ...], ...]:
        """The interaction matrix at a phase's composition.

        The phase-independent rules resolve their matrix once per state and this returns
        it unchanged. **The Soreide-Whitson rule is phase-dependent**: its salinity
        correlation replaces the water-gas entries, and only in a water-rich phase, so the
        matrix is a function of the composition and this is where it is computed.

        The gate is NeqSim's own `water.x > 0.8`: a vapour phase whose water is a few
        molecules in a thousand keeps the base matrix, which is why a model that skipped
        the gate would agree with the oracle on a gas and be wrong on every brine.
        """
        parameters = self.soreide_whitson
        if self.mixing_rule != "soreide_whitson" or parameters is None:
            return self.kij
        roles = parameters.roles
        if not any(
            role == "water" and fraction > 0.8 for role, fraction in zip(roles, x, strict=True)
        ):
            return self.kij
        salinity = parameters.salinity
        return tuple(
            tuple(
                self.kij[i][j]
                if roles[j] != "water"
                else self._aqueous_kij(i, roles[i], reduced.reduced_temperatures[i], salinity)
                for j in range(len(roles))
            )
            for i in range(len(roles))
        )

    def _aqueous_kij(
        self, index: int, role: str, reduced_temperature: float, salinity: float
    ) -> float:
        """One water-gas entry of the Soreide-Whitson correlation.

        `getkijWhitsonSoreideAqueous`'s chain, on the legacy parameterisation - the two
        newer ones (Chabab 2019, Burgoyne-Nielsen 2026) are not ported and this model does
        not claim them.
        """
        if role == "water":
            return 0.0
        if role == "nitrogen":
            return 0.997 * (
                -1.70235 * (1.0 + 0.025587 * math.pow(salinity, 0.75))
                + 0.44338 * (1.0 + 0.08126 * math.pow(salinity, 0.75)) * reduced_temperature
            )
        if role == "carbon_dioxide":
            # NeqSim's ladder has a 0.8 branch above 3.5 mol/kg, but the 0.9 branch above
            # 2.0 fires first, so it is unreachable.
            multip_k = 0.9 if salinity > 2.0 else 1.0
            return (
                multip_k
                * 0.989
                * (
                    -0.31092 * (1.0 + 0.15587 * math.pow(salinity, 0.75))
                    + 0.2358 * (1.0 + 0.17837 * math.pow(salinity, 0.98)) * reduced_temperature
                    - 21.2566 * math.exp(-math.pow(6.7222, reduced_temperature) - salinity)
                )
            )
        omega = self.components[index].omega
        c0, c1, c2 = 0.017407, 0.033516, 0.011478
        a0 = 1.112 - 1.7369 * math.pow(omega, -0.1)
        a1 = 1.1001 + 0.83 * omega
        a2 = -0.15742 - 1.0988 * omega
        return 0.777 * (
            (1.0 + c0 * salinity) * a0
            + (1.0 + c1 * salinity) * a1 * reduced_temperature
            + (1.0 + c2 * salinity) * a2 * reduced_temperature * reduced_temperature
        )


def mixture(
    components: Iterable[Component],
    kij: Mapping[tuple[int, int], float] | None = None,
    cubic: Cubic = PR,
    alpha: str = "pr",
    associating: bool = False,
    mixing_rule: str = "classic",
    umr: UnifacUmrpruParameters | None = None,
    soreide_whitson: SoreideWhitsonParameters | None = None,
    huron_vidal: HuronVidalParameters | None = None,
    furst: FurstElectrolyte | None = None,
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
        associating: whether a phase model runs the Wertheim association
            contribution. See :attr:`Mixture.associating`.
        mixing_rule: ``"classic"``, the default, or ``"umr"``. See
            :attr:`Mixture.mixing_rule`.
        umr: the UNIFAC-UMR-PRU tables the ``"umr"`` rule reads, and ``None`` for
            every other rule. Both of these exist for one model,
            :func:`azoth.eos.components.umr_cpa_mixture_of`, which is the only caller
            that names either.
        soreide_whitson: the roles and salinity the ``"soreide_whitson"`` rule reads,
            and ``None`` for every other rule. The ``kij`` given beside it is that rule's
            base matrix, which is ``INTER.csv``'s ``KIJWhitsonSoriede`` column.

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
        alpha=alpha,
        associating=associating,
        mixing_rule=mixing_rule,
        umr=umr,
        soreide_whitson=soreide_whitson,
        huron_vidal=huron_vidal,
        furst=furst,
    )
