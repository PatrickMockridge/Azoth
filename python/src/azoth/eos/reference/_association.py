"""The Wertheim association contribution, the pure-Python reference.

The twin of ``crates/azoth-eos/src/association.rs``, hand-written like every other
reference here rather than generated from it: the point of two kernels is that they are
compared, and a translation machine would compare a language with itself.

What is here is the **value path** - the site schemes and their bond rule, the
Carnahan-Starling distribution function, the association strength ``Delta``, the
site-fraction solve, and the Helmholtz energy with the fugacity and enthalpy terms. That
is what a *phase state* reads, and it is what ``eos.srk_cpa_phase`` needs.

**What is deliberately not here** is the second-order derivative surface: the implicit
``dV^2``, ``dV dn_j`` and ``dV dT`` families that ``Mixture::phase_derivatives`` reads to
answer a flash's Jacobian. A model that reads ``phase_derivatives`` cannot yet be a CPA
model in either language, and saying so here rather than silently returning the cubic's
derivative is the difference between a missing answer and a wrong one.
:func:`Association.temperature_derivative` is the exception and is first order: the
enthalpy departure is ``A/(RT) - T d(A/RT)/dT``, so the value surface needs it.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

from azoth.core.errors import InvalidInputError, OutOfRangeError

if TYPE_CHECKING:
    from collections.abc import Sequence

#: NeqSim's ``ThermodynamicConstantsInterface.R``.
#:
#: NeqSim's own rounded literal rather than CODATA, for the reason the Rust kernel takes
#: it: the association energy enters ``exp(eps/(R T))``, where a difference in the eighth
#: digit is a difference in the answer.
R = 8.3144621

#: The successive-substitution tolerance, NeqSim's ``solveX2``.
SOLVE_TOLERANCE = 1.0e-12
#: The sweep cap, NeqSim's ``solveX2(15)``.
MAX_SWEEPS = 15


class SiteScheme(Enum):
    """The association site scheme a component carries.

    The names are the values of ``COMP.csv``'s ``associationscheme`` column. A scheme
    fixes both how many sites the molecule has and which of them bond, and the rule is
    NeqSim's ``getInteractionMatrix``: two sites associate exactly when their charges
    differ in sign.
    """

    #: No sites: a component the table gives no scheme, and the stand-in a mixture carries
    #: for every non-associating member so the site offsets line up.
    NON_ASSOCIATING = "none"
    #: One site. Its charge product is positive, so it does not self-associate.
    ONE_A = "1A"
    #: Two sites of the same sign. No pair associates.
    TWO_A = "2A"
    #: Two sites of opposite sign - one donor, one acceptor.
    TWO_B = "2B"
    #: Four sites, two of each sign.
    FOUR_C = "4C"

    @property
    def charges(self) -> tuple[int, ...]:
        """The charges NeqSim gives the sites, ``CPAMixingRuleHandler:32-35``."""
        return _CHARGES[self]

    @property
    def site_count(self) -> int:
        """How many sites the scheme has."""
        return len(_CHARGES[self])

    def associates_within(self, a: int, b: int) -> bool:
        """Whether this scheme's own sites ``a`` and ``b`` associate."""
        charges = _CHARGES[self]
        return charges[a] * charges[b] < 0

    def associates_across(self, a: int, other: SiteScheme, b: int) -> bool:
        """Whether a site on this scheme associates with a site on ``other``."""
        return _CHARGES[self][a] * _CHARGES[other][b] < 0

    @property
    def self_bonds(self) -> bool:
        """Whether this scheme's sites bond with each other at all.

        False for :attr:`ONE_A` and :attr:`TWO_A`, whose charge vectors carry one sign, so
        NeqSim's product test is positive for every pair and the interaction matrix is all
        zeros. Such a component still *cross*-associates with an oppositely-charged
        partner, which is why this is a property of the scheme and not a verdict on the
        component.
        """
        n = self.site_count
        return any(self.associates_within(a, b) for a in range(n) for b in range(n))

    @classmethod
    def from_databank_name(cls, name: str) -> SiteScheme | None:
        """The scheme a databank name denotes, or ``None``.

        ``None`` for ``"0"``, the table's marker for a component with no scheme at all,
        and for any name this library does not carry. A site *count* cannot stand in for
        this: the table carries ``1A`` at zero sites and both ``2A`` and ``2B`` at two, and
        NeqSim's ``setAssociationScheme`` switches on the name.
        """
        return _SCHEME_BY_NAME.get(name.strip())


#: The site charges, NeqSim's `CPAMixingRuleHandler:32-35`. A table keyed by the scheme
#: rather than carried on the enum's value, because the value is the *databank name* and a
#: scheme has two names - the one the table writes and the charges it implies - which are
#: separate things: the table carries `1A` at zero sites.
_CHARGES: dict[SiteScheme, tuple[int, ...]] = {
    SiteScheme.NON_ASSOCIATING: (),
    SiteScheme.ONE_A: (-1,),
    SiteScheme.TWO_A: (-1, -1),
    SiteScheme.TWO_B: (1, -1),
    SiteScheme.FOUR_C: (1, 1, -1, -1),
}

_SCHEME_BY_NAME: dict[str, SiteScheme] = {
    "1A": SiteScheme.ONE_A,
    "2A": SiteScheme.TWO_A,
    "2B": SiteScheme.TWO_B,
    "4C": SiteScheme.FOUR_C,
}


@dataclass(frozen=True, slots=True)
class AssociationComponent:
    """One component's association parameters as the kernel reads them.

    ``energy`` is ``eps`` in J/mol and ``volume`` is ``kappa_AB``, dimensionless. Water's
    16655 J/mol is ``eps/R`` of 2003.1 K against CPA's published 2003.4 K, which is how the
    unit is fixed.
    """

    scheme: SiteScheme
    energy: float
    volume: float


#: A component that carries no association, for a mixture's shape.
NON_ASSOCIATING = AssociationComponent(SiteScheme.NON_ASSOCIATING, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class CrossRule:
    """The cross-association rule for one component pair.

    NeqSim's ``assosSchemeType``: ``0`` is the Elliott rule, which combines two pure
    components' association strengths geometrically, and ``1`` is CR-1, which uses the
    pair's own energy and volume. NeqSim initialises the whole matrix to zero - the Elliott
    rule - and overwrites only the pairs the ``intertemp`` table has a row for.
    """

    #: ``"elliott"`` or ``"cr1"``.
    kind: str
    #: ``eps_ij`` in J/mol, for CR-1.
    energy: float = 0.0
    #: ``kappa_AB`` for the pair, for CR-1.
    volume: float = 0.0


#: The Elliott rule, NeqSim's default for every pair its table has no row for.
ELLIOTT = CrossRule("elliott")

#: Which fitted set a cubic reads, **by geometry and not by name**.
#:
#: `rk` shares Soave's `delta` pair and `pr` Peng-Robinson's, so each reads the family it is
#: shaped like - the same rule `azoth_eos::association::AssociationCubic::of` applies. The
#: two families are separate fits rather than one converted, so this is a selection: water's
#: `kappa_AB` is 0.0692 for SRK against 0.046473789 for PR, and reading the wrong one is a
#: different fluid with no symptom.
_FAMILY_OF_CUBIC: dict[str, str] = {"pr": "pr", "srk": "srk", "rk": "srk"}


def family_of(cubic_name: str) -> str:
    """The associating family a cubic's short name belongs to.

    Raises:
        InvalidInputError: if the name is no cubic this build has.
    """
    try:
        return _FAMILY_OF_CUBIC[cubic_name]
    except KeyError:
        raise InvalidInputError(
            "eos", f"unknown cubic {cubic_name!r}; expected 'pr', 'srk' or 'rk'"
        ) from None


class Rdf:
    """Carnahan-Starling ``g = (1 - eta/2)/(1 - eta)**3`` with ``eta = B/(4V)``.

    ``PhaseCPAInterface.calc_g``, ``calc_lngV`` and ``calc_lngVV``. It is a property of the
    *mixture* - ``g`` carries no component index - so one value serves every component, and
    NeqSim's ``calc_lngi(i)`` is the derivative with respect to ``n_i`` rather than a
    per-component quantity.
    """

    __slots__ = (
        "d2_ln_g_d_eta2",
        "d2_ln_g_dv2",
        "d_ln_g_dn",
        "d_ln_g_dv",
        "eta",
        "g",
    )

    def __init__(self, covolumes: Sequence[float], moles: Sequence[float], v: float) -> None:
        big_b = sum(b * n for b, n in zip(covolumes, moles, strict=True))
        eta = big_b / (4.0 * v)
        one_minus = 1.0 - eta
        self.eta = eta
        self.g = (1.0 - eta / 2.0) / (one_minus * one_minus * one_minus)
        # d ln g / d eta = 3/(1 - eta) - 1/(2 - eta).
        d_ln_g_d_eta = 3.0 / one_minus - 1.0 / (2.0 - eta)
        # d^2 ln g / d eta^2. Both terms are positive: `ln g` carries `-3 ln(1 - eta)`,
        # whose second derivative is `+3/(1 - eta)^2`, and `1 - eta/2` contributes
        # `-1/(2 - eta)` and then `-1/(2 - eta)^2`.
        self.d2_ln_g_d_eta2 = 3.0 / (one_minus * one_minus) - 1.0 / ((2.0 - eta) * (2.0 - eta))
        self.d_ln_g_dn = [b * d_ln_g_d_eta / (4.0 * v) for b in covolumes]
        # eta = B/(4V), so d eta/dV = -B/(4V^2).
        d_eta_dv = -big_b / (4.0 * v * v)
        self.d_ln_g_dv = d_ln_g_d_eta * d_eta_dv
        self.d2_ln_g_dv2 = self.d2_ln_g_d_eta2 * d_eta_dv * d_eta_dv + d_ln_g_d_eta * big_b / (
            2.0 * v * v * v
        )

    @property
    def d_ln_g_d_eta(self) -> float:
        """``d ln g / d eta``, recomputed from the stored packing fraction."""
        return 3.0 / (1.0 - self.eta) - 1.0 / (2.0 - self.eta)


def delta_nog(
    a: AssociationComponent,
    covolume_a: float,
    c: AssociationComponent,
    covolume_c: float,
    rule: CrossRule,
    t: float,
) -> float:
    """The association strength between two components, without the distribution function.

    NeqSim's ``deltaNog``. The covolumes are arguments because they are not association
    parameters - for CPA they are the ``bCPA`` values NeqSim substitutes for the cubic's
    own ``b``, and for PC-SAFT the segment diameter enters differently. Under either rule
    the self-pair ``(i, i)`` reduces to component ``i``'s own strength, so a component's
    self-association needs no special case.
    """

    def strength(x: AssociationComponent, covolume: float) -> float:
        return (math.exp(x.energy / (R * t)) - 1.0) * covolume * x.volume

    if rule.kind == "elliott":
        return math.sqrt(strength(a, covolume_a) * strength(c, covolume_c))
    # CR-1, with the pair's own energy and volume.
    return (math.exp(rule.energy / (R * t)) - 1.0) * (covolume_a + covolume_c) / 2.0 * rule.volume


def _boltzmann_d_ln_dt(energy: float, t: float) -> float:
    """``d ln(exp(eps/(R T)) - 1) / dT``, the only temperature dependence ``Delta`` has."""
    ratio = energy / (R * t)
    boltzmann = math.exp(ratio) - 1.0
    if abs(boltzmann) < 1.0e-50:
        return 0.0
    return -ratio * math.exp(ratio) / (t * boltzmann)


def _delta_nog_d_ln_dt(
    a: AssociationComponent, c: AssociationComponent, rule: CrossRule, t: float
) -> float:
    """The logarithmic temperature derivative of ``DeltaNog`` for a pair."""
    if rule.kind == "elliott":
        # `sqrt(K_a K_c)`, so the log-derivative is the mean of the two components'.
        return 0.5 * (_boltzmann_d_ln_dt(a.energy, t) + _boltzmann_d_ln_dt(c.energy, t))
    return _boltzmann_d_ln_dt(rule.energy, t)


@dataclass(frozen=True, slots=True)
class SiteState:
    """The converged site fractions and the energies they carry."""

    #: ``X_A``, one entry per site, flattened in component order.
    fractions: tuple[float, ...]
    #: Successive-substitution sweeps taken, including the one that met the tolerance.
    iterations: int
    #: Whether the solve met :data:`SOLVE_TOLERANCE` before :data:`MAX_SWEEPS`.
    converged: bool
    #: Whether the Newton refinement ran and left every fraction positive.
    refined: bool
    #: ``hCPA``, the unbonded site count weighted by moles: ``sum_i n_i sum_{A in i}(1 - X_A)``.
    #:
    #: NeqSim's ``PhaseCPAInterface.calc_hCPA``, and the factor multiplying the distribution
    #: function's derivative in every association term. It is *not* one, which is what
    #: comparing against the probe settled: the port first carried the class field's
    #: initialiser and was 4% out on ``dFCPAdN``.
    unbonded_sites: float
    #: ``A_assoc / (R T)``, the association Helmholtz energy over ``R T``.
    helmholtz_rt: float
    #: ``d (A_assoc/RT) / dV`` at constant composition and temperature.
    d_helmholtz_dv: float
    #: The association contribution to ``ln phi_i``, one entry per component, at constant
    #: ``T`` and ``V``. NeqSim's ``ComponentSrkCPA.dFCPAdN``.
    ln_phi: tuple[float, ...]


class Association:
    """The whole mixture's association parameters."""

    __slots__ = ("components", "cross", "offsets")

    def __init__(
        self, components: Sequence[AssociationComponent], cross: Sequence[CrossRule] = ()
    ) -> None:
        """A mixture's association parameters.

        Raises:
            InvalidInputError: if ``cross`` is neither empty nor ``N * N`` long.
            OutOfRangeError: if any component's energy or volume is negative. Both enter an
                exponential and a square root, and a negative volume would make the Elliott
                combination the square root of a negative number.
        """
        n = len(components)
        if cross and len(cross) != n * n:
            raise InvalidInputError(
                "cross",
                f"a cross rule matrix for {n} components has {len(cross)} entries",
            )
        for i, component in enumerate(components):
            for field_name, value in (
                ("association energy", component.energy),
                ("association volume", component.volume),
            ):
                if not math.isfinite(value) or value < 0.0:
                    raise OutOfRangeError(
                        field_name,
                        value,
                        f"component {i}'s {field_name} is non-negative by definition: it "
                        f"enters an exponential and a square root, so a negative one is not "
                        f"a weaker association but no association at all",
                    )
        self.components = tuple(components)
        self.cross = tuple(cross) if cross else (ELLIOTT,) * (n * n)
        offsets = [0]
        total = 0
        for component in self.components:
            total += component.scheme.site_count
            offsets.append(total)
        self.offsets = tuple(offsets)

    def __len__(self) -> int:
        """The number of components."""
        return len(self.components)

    @property
    def site_count(self) -> int:
        """The total number of association sites over every component."""
        return self.offsets[len(self.components)]

    def sites_of(self, i: int) -> range:
        """The sites belonging to component ``i``, as a range into the flattened vector."""
        return range(self.offsets[i], self.offsets[i + 1])

    def component_of_site(self, site: int) -> int:
        """Which component a flattened site index belongs to."""
        return bisect_right(self.offsets, site) - 1

    def cross_rule(self, i: int, j: int) -> CrossRule:
        """The cross rule for a component pair."""
        return self.cross[i * len(self.components) + j]

    def bonds(self, a: int, b: int) -> bool:
        """Whether the sites ``a`` and ``b`` bond, across components if they differ.

        The pair's :class:`CrossRule` is not consulted: it changes a bond's strength, never
        whether one exists.
        """
        ia, ib = self.component_of_site(a), self.component_of_site(b)
        sa, sb = a - self.offsets[ia], b - self.offsets[ib]
        if ia == ib:
            return self.components[ia].scheme.associates_within(sa, sb)
        return self.components[ia].scheme.associates_across(sa, self.components[ib].scheme, sb)

    def has_bonds(self) -> bool:
        """Whether any pair of sites in the mixture bonds at all.

        A mixture of ``ONE_A`` and ``TWO_A`` components has sites and no bonds, so every
        site fraction is one and the association energy is exactly zero.
        """
        n = self.site_count
        return any(self.bonds(a, b) for a in range(n) for b in range(n))

    def delta_matrix(self, covolumes: Sequence[float], t: float, rdf: Rdf) -> list[float]:
        """``Delta_ij = bond(i,j) deltaNog_ij g``, row-major.

        NeqSim builds this once per state in ``initCPAMatrix`` and reuses it for the solve
        and for every derivative, which is why it is returned rather than recomputed.
        """
        n = self.site_count
        delta = [0.0] * (n * n)
        for i in range(n):
            ci = self.component_of_site(i)
            for j in range(n):
                if not self.bonds(i, j):
                    continue
                cj = self.component_of_site(j)
                delta[i * n + j] = (
                    delta_nog(
                        self.components[ci],
                        covolumes[ci],
                        self.components[cj],
                        covolumes[cj],
                        self.cross_rule(ci, cj),
                        t,
                    )
                    * rdf.g
                )
        return delta

    def solve(
        self, covolumes: Sequence[float], moles: Sequence[float], v: float, t: float
    ) -> SiteState:
        """Solve for the site fractions and evaluate the association energy.

        NeqSim's ``solveX``: successive substitution to :data:`SOLVE_TOLERANCE`, then a
        Newton refinement on the same residual. The two are not redundant - the
        substitution is globally reliable and converges linearly, the Newton is quadratic
        and can leave the physical branch - so both are kept and each reports whether it
        finished.

        Raises:
            InvalidInputError: if the state vectors are not one entry per component.
            OutOfRangeError: if ``t`` or ``v`` is not positive.
        """
        n = len(self.components)
        for field_name, given in (("covolumes", len(covolumes)), ("moles", len(moles))):
            if given != n:
                raise InvalidInputError(
                    field_name, f"a mixture of {n} components takes {given} entries"
                )
        for field_name, value in (("v", v), ("T", t)):
            if not math.isfinite(value) or value <= 0.0:
                raise OutOfRangeError(
                    field_name,
                    value,
                    "the association solve needs a positive volume and temperature",
                )

        sites = self.site_count
        if sites == 0:
            return SiteState(
                fractions=(),
                iterations=0,
                converged=True,
                refined=False,
                unbonded_sites=0.0,
                helmholtz_rt=0.0,
                d_helmholtz_dv=0.0,
                ln_phi=(0.0,) * n,
            )

        rdf = Rdf(covolumes, moles, v)
        delta = self.delta_matrix(covolumes, t, rdf)
        # `K_ij = n_j Delta_ij / V`, so that `sum_j K_ij X_j` is `S_i` and the fixed point
        # `X_i = 1/(1 + S_i)` is NeqSim's `solveX2`.
        #
        # NeqSim's *Newton* path builds `Klk_ij = n_i n_j Delta_ij / V` instead, but pairs
        # it with the `n_i`-scaled residual `n_i (1/X_i - 1) - sum_j Klk_ij X_j`, which
        # vanishes at the same point. Taking that matrix without its scaling gives each site
        # fraction a spurious factor of its component's mole number, which is what a finite
        # difference against the substitution catches and nothing else does.
        klk = [
            moles[self.component_of_site(j)] * delta[i * sites + j] / v
            for i in range(sites)
            for j in range(sites)
        ]

        # Successive substitution, Gauss-Seidel: the site written this sweep is visible to
        # the rest of the sweep, which is NeqSim's in-place update.
        x = [1.0] * sites
        iterations = 0
        converged = False
        while iterations < MAX_SWEEPS:
            iterations += 1
            error = 0.0
            for i in range(sites):
                old = x[i]
                new = 1.0 / (1.0 + sum(klk[i * sites + j] * x[j] for j in range(sites)))
                x[i] = new
                error += abs((old - new) / new)
            if error < SOLVE_TOLERANCE:
                converged = True
                break

        refined = _newton_refine(x, klk, sites)

        # `A_assoc/(RT) = sum_i n_i sum_{A in i} (ln X_A - X_A/2 + 1/2)`.
        helmholtz_rt = sum(
            n_i * sum(math.log(x[a]) - x[a] / 2.0 + 0.5 for a in self.sites_of(i))
            for i, n_i in enumerate(moles)
        )
        # `hCPA = sum_i n_i sum_{A in i} (1 - X_A)`: the number of unbonded sites, weighted
        # by the moles that carry them.
        unbonded_sites = sum(
            n_i * sum(1.0 - x[a] for a in self.sites_of(i)) for i, n_i in enumerate(moles)
        )
        ln_phi = tuple(
            sum(math.log(x[a]) for a in self.sites_of(i)) - 0.5 * unbonded_sites * rdf.d_ln_g_dn[i]
            for i in range(n)
        )
        # NeqSim's `dFCPAdV = (hCPA/(2V)) (1 - V d ln g/dV)`.
        d_helmholtz_dv = unbonded_sites * (1.0 - v * rdf.d_ln_g_dv) / (2.0 * v)

        return SiteState(
            fractions=tuple(x),
            iterations=iterations,
            converged=converged,
            refined=refined,
            unbonded_sites=unbonded_sites,
            helmholtz_rt=helmholtz_rt,
            d_helmholtz_dv=d_helmholtz_dv,
            ln_phi=ln_phi,
        )

    def temperature_derivative(
        self,
        covolumes: Sequence[float],
        moles: Sequence[float],
        v: float,
        t: float,
        state: SiteState,
    ) -> float:
        """``d(A_assoc/(R T))/dT`` at constant composition and volume.

        The one *first*-order implicit derivative the value surface needs, and the only one
        here: the enthalpy departure is ``A/(RT) - T d(A/RT)/dT``, so a model that reports
        ``h_dep`` cannot do without it. The site fractions solve a nonlinear system, so
        this is one linear solve against that system's Jacobian rather than a formula.

        The residual is ``F_i = X_i (1 + S_i) - 1`` with ``S_i = (1/V) sum_k m_k Delta_ik X_k``,
        and ``T`` moves ``DeltaNog`` alone - the distribution function carries no
        temperature. So the right-hand side is ``X_i`` times the ``T``-derivative of ``S_i``.
        """
        sites = self.site_count
        if sites == 0:
            return 0.0
        rdf = Rdf(covolumes, moles, v)
        delta = self.delta_matrix(covolumes, t, rdf)
        x = state.fractions

        s = [
            sum(
                moles[self.component_of_site(k)] * delta[i * sites + k] * x[k] for k in range(sites)
            )
            / v
            for i in range(sites)
        ]
        jacobian = [
            (1.0 if i == k else 0.0) * (1.0 + s[i])
            + x[i] * moles[self.component_of_site(k)] * delta[i * sites + k] / v
            for i in range(sites)
            for k in range(sites)
        ]
        # Negated, because the step solves `J X' = -dF/dT`.
        rhs = [
            -x[i]
            * sum(
                moles[self.component_of_site(k)]
                * delta[i * sites + k]
                * _delta_nog_d_ln_dt(
                    self.components[self.component_of_site(i)],
                    self.components[self.component_of_site(k)],
                    self.cross_rule(self.component_of_site(i), self.component_of_site(k)),
                    t,
                )
                * x[k]
                for k in range(sites)
            )
            / v
            for i in range(sites)
        ]
        if not _solve_many(jacobian, [rhs], sites):
            raise InvalidInputError(
                "site fractions",
                "the site-fraction Jacobian is singular at this state, so the association's "
                "temperature derivative is not defined here",
            )
        return sum(
            moles[self.component_of_site(a)] * (1.0 / x[a] - 0.5) * rhs[a] for a in range(sites)
        )


def _newton_refine(x: list[float], klk: Sequence[float], sites: int) -> bool:
    """Newton-refine the site fractions in place, returning whether it stayed physical.

    NeqSim's ``solveX`` after ``solveX2``. The residual is
    ``F_i = X_i (1 + sum_j Klk_ij X_j) - 1``, whose fixed point is the substitution, and the
    step is ``-J^{-1} F``. A fraction that would go non-positive is clipped to ``1e-10`` -
    NeqSim's value - which leaves the solve unrefined rather than wrong, so the caller keeps
    the substitution's answer.
    """
    if sites == 0:
        return False
    newton_tolerance = 1.0e-12
    max_newton = 100
    physical = True
    for _ in range(max_newton):
        residual = [0.0] * sites
        jacobian = [0.0] * (sites * sites)
        for i in range(sites):
            inner = sum(klk[i * sites + j] * x[j] for j in range(sites))
            residual[i] = x[i] * (1.0 + inner) - 1.0
            for k in range(sites):
                jacobian[i * sites + k] = (1.0 if i == k else 0.0) * (1.0 + inner) + x[i] * klk[
                    i * sites + k
                ]
        if all(abs(r) < newton_tolerance for r in residual):
            return physical
        # One right-hand side, so a fresh factorization each step - the Jacobian moves with
        # `X` and reusing it would solve the previous step's system.
        step = [list(residual)]
        if not _solve_many(jacobian, step, sites):
            return False
        for i in range(sites):
            x[i] -= step[0][i]
            if x[i] <= 0.0:
                x[i] = 1.0e-10
                physical = False
    return physical


def _solve_many(a: Sequence[float], rhs: list[list[float]], n: int) -> bool:
    """Solve ``A x = b`` for several right-hand sides at once, in place.

    The association kernel's Jacobian is the only dense solve in this reference and it is
    ``sites x sites`` - a handful - which is why it is written here rather than taken as a
    dependency. The factorization is shared across the right-hand sides because the implicit
    derivatives are one solve each against the *same* matrix, which is also why ``a`` is not
    consumed: elimination runs on a copy.

    Returns ``False`` on a singular matrix, leaving ``rhs`` untouched, so a caller reports an
    undefined derivative rather than a wrong one.
    """
    lu = list(a)
    for col in range(n):
        pivot = col
        for row in range(col + 1, n):
            if abs(lu[row * n + col]) > abs(lu[pivot * n + col]):
                pivot = row
        if lu[pivot * n + col] == 0.0:
            return False
        if pivot != col:
            for k in range(n):
                lu[col * n + k], lu[pivot * n + k] = lu[pivot * n + k], lu[col * n + k]
            for column in rhs:
                column[col], column[pivot] = column[pivot], column[col]
        for row in range(col + 1, n):
            factor = lu[row * n + col] / lu[col * n + col]
            if factor == 0.0:
                continue
            for k in range(col, n):
                lu[row * n + k] -= factor * lu[col * n + k]
            for column in rhs:
                column[row] -= factor * column[col]
    for column in rhs:
        for row in range(n - 1, -1, -1):
            total = column[row]
            for k in range(row + 1, n):
                total -= lu[row * n + k] * column[k]
            column[row] = total / lu[row * n + row]
    return True


__all__ = [
    "ELLIOTT",
    "MAX_SWEEPS",
    "NON_ASSOCIATING",
    "SOLVE_TOLERANCE",
    "Association",
    "AssociationComponent",
    "CrossRule",
    "R",
    "Rdf",
    "SiteScheme",
    "SiteState",
    "delta_nog",
]
