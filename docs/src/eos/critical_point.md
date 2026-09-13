# Mixture critical point

`eos.critical_point`

The temperature, pressure and volume at which a mixture of fixed composition stops having two phases to distinguish, from the two criticality conditions of Heidemann and Khalil (1980). A nested Newton: an inner solve drives the smallest eigenvalue of the scaled Helmholtz Hessian to zero at a fixed volume, and an outer solve drives the cubic form along the eigenvector of that eigenvalue to zero.
A *procedure* rather than a direct model: there is an iteration here, a starting point, and two tolerances, and every one of them is part of what is being claimed.

## Source

**Heidemann, R. A.; Khalil, A. M. (1980)** (The calculation of critical points. AIChE Journal 26(5), 769-779. The paper is the method's source and has **not been read**. The implementation follows an open-source one instead, and `verification.notes` records where the two differ.
) - TODO: source needed

DOI: [TODO: source needed](https://doi.org/TODO: source needed)

**Unverified.** The equation is standard, but its citation has not been checked against the primary source by a person.

The method is Heidemann and Khalil's; what is unconfirmed is the equation number and the reading of the paper.
# Where this differs from NeqSim's implementation
The two conditions are the same, and the nested structure is the same. Four things differ, and each was a decision rather than a transcription.
**The derivative frame is the Helmholtz one directly.** NeqSim forms `Q` from `getdfugdn`, which is a constant-pressure composition derivative, and converts to constant volume with `-dfugdp * V_j * dP/dV`. This implementation never forms the constant-pressure derivative: `_mixture_state.criticality_matrix` builds the Hessian of the Helmholtz energy at constant temperature and volume, so there is nothing to convert. That also makes the symmetry structural, and NeqSim needs an explicit symmetrisation step because theirs is not.
**The pressure comes from the volume, not from the cubic.** NeqSim sets the volume and asks its phase object for the pressure through the cubic solver. A critical point is exactly where a cubic is degenerate - at `Tr = Pr = 1` the three roots have merged and even the full-precision constants leave about five significant digits - so this implementation evaluates the explicit equation of state instead, `P = R T/(V - b) - a/(V^2 + 2 b V - b^2)`, which is well posed everywhere the state is.
**The convergence target is the eigenvalue, not a determinant.** NeqSim's own code comment records that a determinant is "far better scaled" as a Rayleigh quotient, and its implementation does use the quotient; it keeps the name `detM` for the variable and reports it as such. This implementation names the quantity for what it is, and takes the **algebraically smallest** eigenvalue rather than the smallest in magnitude, because the one that crosses zero at a critical point is the one that goes from negative to positive.
**The ideal part of the cubic form is written down.** NeqSim's `calcdpd` forms a central second difference of the Rayleigh quotient, which measures the quartic term; the criticality condition is on the cubic one, which is the central *first* difference. This implementation takes the first difference and subtracts the ideal part's third derivative explicitly, without which a pure component reads `1.000000000` at its own critical point instead of zero.
**No validation was inherited, because there is none to inherit.** The NeqSim class checks its result against nothing. Every check below is this project's.
The two conditions are not ours: that a critical point is where the smallest eigenvalue of the scaled Helmholtz Hessian vanishes, and where the cubic form along its eigenvector vanishes, is Heidemann and Khalil's result. What is ours and unverified is the derivative frame, the pressure route, the initialisation, the two tolerances, and the finite-difference step the cubic form is evaluated at.
**The paper has not been read**, so this ships `unverified` even though it is checked hard against closed-form answers. Reading Heidemann and Khalil is what would move it, and it is the same human-blocked work item the README records.
# What it is checked against, and it is not another implementation
A pure component's critical point is known **in closed form** for Peng-Robinson: `Tr = Pr = 1` with `Z_c = (1 - omega_b)/3 = 0.30740130869870386`. So the whole construction - the Hessian at constant volume, the ideal part, the scaling, the cubic form, the nesting - is checked against an analytic answer rather than against a second implementation of the same reading.
Measured on **both** implementations: `Tc` to 2.5e-12 relative, `Pc` to 1.6e-11, and `Z_c` to 5.2e-10 absolute - the same four figures for propane, methane, n-butane and carbon dioxide, which is the signature of a systematic floor rather than of noise. The floor is the finite-difference step the cubic form is evaluated over, and it is what the case tolerance of `1e-9` is set from.
# The test that would catch a plausible wrong answer
The mechanical route - solving `dP/dV = d2P/dV2 = 0` at fixed composition - is exact for a pure component and reproduces that same `0.3074`. It is exact because in reduced variables it has a **single universal root**, so it returns the same `Z_c` for every mixture there is. A pure-component check cannot tell the two routes apart, and a model that only had one would pass.
What tells them apart is that **a mixture's `Z_c` varies with composition**. Measured here for methane / n-butane at `kij = 0.05`: `0.39696`, `0.47865`, `0.54247` and `0.49913` at methane fractions 0.2, 0.4, 0.6 and 0.8 - a spread of 0.146 against the mechanical route's constant. The `Tc` locus falls monotonically from 411.1 K to 273.5 K between the pure endpoints of 425.12 K and 190.56 K, which is the shape a binary's critical locus has.
# Where this is least trustworthy
**The cubic form is a finite difference**, and it is the weaker of the two conditions. Its step is fixed at `1e-4` in the composition, which puts a floor of roughly `1e-10` on the residual that no amount of iterating improves - so the outer tolerance is `1e-9` and the reported residual sits near that floor rather than at it. An analytic third derivative would remove the floor; it would also be a much larger expression, and the finite difference is checked to be converging by the pure-component result rather than assumed to be.
**`Vc` is reported from the iteration, not from the cubic.** It does not come back as `Z_c R Tc/Pc` to more than the iteration's tolerance, and the two agree to about 1e-11 relative for a pure component.

## Algorithm

A model rather than a calculation: what this page pins down is the procedure,
not an equation, and both implementations read it from here.

| Setting | Value |
|---|---|
| Scheme | `heidemann_khalil_critical` |
| Convergence | `absolute` |
| Tolerance | `1e-09` |
| Max iterations | `60` |
| Initialisation | `kay_rule_and_covolume` |

## Inner procedure

| Setting | Value |
|---|---|
| Scheme | `newton_on_temperature_eigenvalue` |
| Convergence | `absolute` |
| Tolerance | `1e-10` |
| Max iterations | `60` |

## Inputs

| Name | Unit | Description |
|---|---|---|
| `Tc` | K | critical temperatures, in the mixture's component order |
| `Pc` | Pa | critical pressures, in the same order |
| `omega` | dimensionless | acentric factors, in the same order. All three vectors are the caller's - this library ships no component databank - and they must be mutually consistent, which is not checked. |
| `kij` | dimensionless | binary interaction parameters, as in `eos.pt_flash` |
| `z` | dimensionless | the composition whose critical point is wanted, checked rather than renormalised. Unlike the flash's `z` this is not a feed being split: it is the composition of the single phase that is about to stop existing. |


## Outputs

| Name | Unit | Description |
|---|---|---|
| `Tc` | K | the critical temperature |
| `Pc` | Pa | the critical pressure, evaluated from the equation of state at the returned volume |
| `Vc` | m**3/mol | the critical molar volume |
| `Z_c` | dimensionless | `Pc Vc/(R Tc)`. Reported because it is the quantity that distinguishes this model from the mechanical route: a pure component's is `(1 - omega_b)/3` and a mixture's varies with composition, which the mechanical conditions cannot produce. |

| Bound | On violation | Why |
|---|---|---|
| `Tc > 0` | raises | the critical temperature of the composition searched for. The iteration is perfectly capable of leaving the region where the equation of state has a critical point and converging on a temperature that is not a state, and this is the check that says so rather than returning it. |
| `Pc > 0` | raises | the pressure at the returned volume. It is evaluated from the explicit equation of state rather than solved for, so a non-positive value means the volume is inside the co-volume and the state is not one. |
| `Z_c > 0` | raises | `Pc Vc/(R Tc)`. **Not an accuracy bound and not a physical range check** - a cubic's critical compressibility is the equation's, not the substance's. It is here because a sign error in any of the three quantities above shows up as a negative `Z_c` before it shows up anywhere else. |

## Assumptions

- **`Tc`, `Pc`, `omega` and `kij` are the caller's and their correctness is NOT CHECKED.** This library ships no component databank and no binary-interaction table.
- **the composition has exactly one critical point.** A mixture can have more than one, and this model finds the one nearest its starting point without any check that there is not another. It is the same limitation `eos.pt_flash` records for a second liquid phase, and it has the same cause: no stability analysis.
- the critical point found is a *mixture* critical point of the given composition, not a point on the phase envelope of a reservoir fluid. The two coincide for a binary and need not for a fluid with more components.
- Peng-Robinson's `Z_c` is not the experimental one for any substance. The model reports what the equation of state gives, and the equation of state is the caller's choice.
- the result is a state of one phase. Nothing checks that the composition is one that could exist at that state - a composition inside a two-phase region has a critical point all the same, and this will find it.

## Cases

| Case | Inputs | Expected |
|---|---|---|
| `pure_propane` | Tc = [369.83], Pc = [4248000.0], omega = [0.1523], kij = [[0.0]], z = [1.0] | Tc = 369.83000000000004, Pc = 4248000.000000011, Vc = 0.00022251409546288967, Z_c = 0.3074013091160614 |
| `methane_and_butane` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], z = [0.4, 0.6] | Tc = 389.60393756196225, Pc = 8496149.7632273, Vc = 0.00018249737998221664, Z_c = 0.4786535349108495 |

## References

- Heidemann, R. A.; Khalil, A. M. (1980). "The calculation of critical points." AIChE Journal 26(5), 769-779. DOI 10.1002/aic.690260510. (the two conditions; the citation is unconfirmed and the paper has NOT been read)
- Michelsen, M. L.; Mollerup, J. M. (2004). "Thermodynamic Models: Fundamentals and Computational Aspects", ch. 5. (states the same conditions and the derivation behind them; not read, and the chapter reference is the one NeqSim's documentation cites)
- Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State." Industrial & Engineering Chemistry Fundamentals 15(1), 59-64. DOI 10.1021/i160057a011. (the equation of state the Helmholtz energy is for; the citation is confirmed and the equation numbers are not)
- NeqSim 3.20.0, `CriticalPointFlash.java`, Apache-2.0, commit dedba8735d030c6411e09b6fd7e69f6c4a136114. (the implementation consulted for the method; it validates its result nowhere, so it is a source for the algorithm and not for the answer)
