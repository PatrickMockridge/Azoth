# Pressure-temperature flash

`eos.pt_flash`

The two-phase split of a mixture at a fixed temperature and pressure: how much of it is vapour, what the two phases are made of, and the fugacity coefficients that make it so. Solved by successive substitution from Wilson K-value estimates, with Rachford-Rice solved by bisection each iteration.
This is the model layer's reason for existing. Everything below it in `eos` is a scalar kernel with a worked example a reader can retrace on a calculator; a flash is an iteration whose answer is a composition vector, and what has to be pinned down is not an equation - `K_i = phi_i^L / phi_i^V` is one line - but the procedure that searches for the state where it holds.

## Source

**Rachford, H. H.; Rice, J. D. (1952); Michelsen, M. L. (1982)** (The Rachford-Rice equation is Rachford & Rice, Journal of Petroleum Technology 4(10), 19-3, DOI 10.2118/952327-G, the same paper `eos.rachford_rice_binary` cites. Successive substitution with a fugacity-coefficient update is the standard method and is set out in Michelsen, M. L. "The isothermal flash problem. Part I. Stability." Fluid Phase Equilibria 9(1), 1-19 (1982), DOI 10.1016/0378-3812(82)85001-2.
) - TODO: source needed

DOI: [10.1016/0378-3812(82)85001-2](https://doi.org/10.1016/0378-3812(82)85001-2)

**Unverified.** The equation is standard, but its citation has not been checked against the primary source by a person.

The equations are not in doubt. `K_i = phi_i^L / phi_i^V`, `sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)) = 0`, and the material balance `z = (1 - beta) x + beta y` are the definition of an isothermal two-phase flash, and the second is the equation `eos.rachford_rice_binary` already registers.
What is ours, and unverified, is the **procedure** - the initialisation, the bracket, the tolerances, the stopping rule. No source states our choices because no source would. They are recorded here so a reader judges the choices rather than reverse-engineering them, and the three subsections below are the parts of the procedure that were **wrong in an earlier draft of this spec and are recorded with the measurement that caught them**. That is the content worth reading; the rest is bookkeeping.
# Correction 1: the Rachford-Rice bracket, and what "no root" means
This is the one that took three attempts, and the first two were wrong in the same direction - each was a *patch* on the textbook statement rather than a derivation from it.
The textbook bracket for `g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))` is `1/(1 - K_max) < beta < 1/(1 - K_min)`. That is correct **only when the K-values straddle one**. `g` has poles at `1/(1 - K_i)`, is strictly decreasing between consecutive poles, and jumps from `-inf` to `+inf` across each. That fixes the whole root structure:
* the leftmost interval runs from `0-` at `-inf` to `-inf` at the first pole -
  **no root**;
* the rightmost runs from `+inf` at the last pole to `0+` at `+inf` -
  **no root**;
* every bounded interval between consecutive poles has **exactly one** root, and
  there are `n - 1` of them.

Of those `n - 1` roots, exactly one has both phases positive, and it is the one on which `1 + beta (K_i - 1) > 0` for every `i` - because that is what keeps `x_i = z_i / (1 + beta (K_i - 1))` in `[0, 1]`. Solving those inequalities gives the interval the bisection must search:
```text lower = max over {i : K_i > 1} of 1/(1 - K_i) upper = min over {i : K_i < 1} of 1/(1 - K_i) ```
and **it exists only when the K-values straddle one.** When every `K_i` is on the same side of 1 there is no such interval, because the two expressions are then both negative or both positive - and this is not a numerical failure but a proof. `sum_i y_i = sum_i K_i x_i = 1` with every `K_i < 1` and `sum_i x_i = 1` requires `1 < 1`; with every `K_i > 1` it requires `1 > 1`. Either way the feed has no two-phase solution at those K-values and is single phase.
# What the two wrong versions did
The first draft used the textbook bracket unguarded. Measured: an all-liquid feed (methane/n-butane, T = 200 K, P = 50 MPa, z = (0.6, 0.4)) produced a *negative molar composition* - `b_reduced = -4.78`, which `eos.pr_z_factor` refused - and an all-vapour feed at 1 bar produced `beta = -0.129` from a bracket of `(-0.0022, 0)`.
The second draft replaced it with `min(0, 1/(1 - K_max))` and `max(0, 1/(1 - K_min))`, which is right whenever the K-values straddle one and wrong whenever they do not: it substitutes the interval *containing zero* for the interval *containing the root*, and those differ exactly when every K is on one side. Measured, with propane/n-butane at T = 280 K, P = 1 bar, z = (0.1, 0.9): `K = (5.92, 1.29)`, both above one, so there is no root at all - and the bisection walked to the bracket's end and returned `beta = -2.9e-15` with a vapour mole fraction of **1.162**. A composition above one is the kind of answer a caller can only catch by checking, and nothing in the model was checking.
Both are now tested. `both_compositions_stay_positive` in both languages sweeps two binaries over 60 states and asserts every returned mole fraction is in `[0, 1]`, which is the property the correct bracket is *equivalent* to rather than merely consistent with.
# Correction 2: what the trivial solution is
The plan for this milestone defined `TRIVIAL_SOLUTION` as convergence to `x = y = z` "while the RR solution lay strictly inside (0, 1)". The second clause is wrong, and measurably so. At the trivial solution `g` is identically zero - every `K_i` is 1 - so **`beta` is indeterminate**: it is not inside (0, 1), not outside it, but undefined. Successive substitution approaches the point through geometrically growing `beta`, and where the bisection stops is a ratio of round-off.
Measured: methane/n-butane at T = 330 K, P = 1 bar, z = (0.6, 0.4) - a feed that is unambiguously all *vapour* - converged with `max|ln K_i| = 8.7e-19` and `beta = -7.7e10`. A guard written on `beta < 0` would have labelled it "all liquid". The guard is on the K-values alone:
```text max_i |ln K_i| < 1e-08   during or after convergence  =>  trivial ```
and `beta` is then reported as **absent**, not as a number. That is the one place this model's result shape differs from every other result in the library, and it is deliberate: a fabricated `beta` would be indistinguishable from a real one, and the whole point is that a caller must not be able to mistake this outcome for a phase split.
That particular feed is no longer a trivial solution - with the bracket correction 1 describes, its K-values are all above one and it is *proved* single phase vapour on the first iteration, with `beta` absent for a better reason. The measurement stands as the reason the guard is not written on `beta`, and the `trivial` outcome is still reached where the K-values straddle one throughout: the same mixture and composition at 20 MPa, in 59 iterations, where the feed is a compressed liquid the model can only say it cannot split.
# Correction 3: one of the plan's invariants is vacuous
The plan's §2.5 listed `ln K_i = ln phi_i^L - ln phi_i^V` as an invariant to test to 1e-10 on every case. For a successive-substitution flash that is the **update rule**, not a check on it - the iteration assigns exactly that - so the test can only ever measure the convergence residual, which the result already reports. It is kept in the test suite for what it does say (the loop stopped where it claimed), but the invariants that actually discriminate a right answer from a plausible wrong one are the material balance, the normalisation, `K_i = y_i / x_i`, the Rachford-Rice residual, and the Gibbs minimum. Those five are the ones worth sabotaging.
# Correction 4: `residual` is not a value a case can pin
The first draft of the cases below asserted `residual` alongside the rest, at the case's 1e-9 relative tolerance. Two implementations agreeing to the last bit on `beta`, `x`, `y` and `K` disagreed on it in the sixth significant figure - 4.658240e-11 against 4.658231e-11.
That is not a defect. `residual` is an rms difference of quantities built from `ln` and `sqrt`, neither of which is correctly rounded, so two implementations differ in its last few *ulps*, and a value of 4.7e-11 has no significant figures left to spare. The relative tolerance the case declares is meaningless at that magnitude, and pinning the number would have been a claim about the C library rather than about this model.
So the cases pin what is reproducible - `beta`, both compositions, the K-values, both `ln phi` vectors, both roots, `min_t_over_tc`, and the **iteration count**, which is the sharp structural check and is exactly reproducible - and `residual` is asserted against the algorithm's own tolerance instead, in the test files. That is the claim that matters about it, and the one that is true.
# The honest gap: there is no stability test
Successive substitution finds *a* stationary point of the flash equations. It does not ask whether the feed was stable in the first place, and the two are different questions. Correction 1 closes part of the gap - a feed with no Rachford-Rice root is now *proved* single phase rather than approximated as one, and that is a diagnostic the textbook bracket cannot make - but it is not a stability test and the difference is measured.
What remains is the `trivial` outcome, which is the gap stated plainly: the iteration converges to `x = y = z`, the feed is single phase, and the model **cannot say which phase**. For methane/n-butane at z = (0.6, 0.4) the trivia outcome is reached at 20 MPa across 280-330 K, where the feed is a compressed liquid, and the model knows that only in the sense that it cannot find a split. A tangent-plane analysis is what answers it, and it is a later milestone. A caller who needs to know must not read `phase: trivial` as "all liquid".
The three states that *are* diagnosed, for the same mixture at the same composition, over a pressure sweep at 330 K:
| P | phase | beta | iterations | |---|---|---|---| | 100 kPa | `all_vapour` | absent | 1 | | 1 MPa | `all_vapour` | 1.674 | 9 | | 2.5-10 MPa | `two_phase` | 0.845 to 0.470 | 13-42 | | 20 MPa | `trivial` | absent | 59 | | 50 MPa | `all_liquid` | absent | 1 |
Three distinct routes to a single phase - no root, a negative flash, and the trivial convergence - and only the third is a failure to know.
# No accuracy claim
Peng-Robinson's VLE predictions are good to a few per cent for light hydrocarbons and worse for polar and asymmetric mixtures. As in the rest of this namespace, no bound below is an accuracy claim - they are all about where the calculation is defined.

## Algorithm

A model rather than a calculation: what this page pins down is the procedure,
not an equation, and both implementations read it from here.

| Setting | Value |
|---|---|
| Scheme | `successive_substitution_flash` |
| Convergence | `absolute` |
| Tolerance | `1e-10` |
| Max iterations | `300` |
| Initialisation | `wilson` |

## Inner procedure

| Setting | Value |
|---|---|
| Scheme | `rachford_rice_bisection` |
| Convergence | `absolute` |
| Tolerance | `1e-14` |
| Max iterations | `200` |

## Inputs

| Name | Unit | Description |
|---|---|---|
| `Tc` | K | critical temperatures, in the mixture's component order |
| `Pc` | Pa | critical pressures, in the same order |
| `omega` | dimensionless | acentric factors, in the same order. All three of these vectors are the caller's, though `azoth.eos.component` will look them up - and they must be mutually consistent, which is not checked. |
| `kij` | dimensionless | binary interaction parameters. Zero diagonal because a component does not interact with itself, and the implementation reads only the upper triangle - a caller who supplies an asymmetric matrix is not corrected. |
| `T` | K | absolute temperature |
| `P` | Pa | absolute pressure |
| `z` | dimensionless | overall mole fractions. Checked rather than renormalised: the entries must be non-negative and sum to one to within a tolerance, and a feed that does not is refused. Silently renormalising would make a caller's composition error invisible in a way that changes every number downstream. |


## Outputs

| Name | Unit | Description |
|---|---|---|
| `beta` | dimensionless | *Optional.* The vapour fraction, in `[0, 1]` for a two-phase split. Outside `[0, 1]` the feed is single phase and this is the negative-flash value - the amount of the absent phase that would have to be added to bring the feed to saturation - reported with `OUT_OF_VALID_RANGE` rather than suppressed, because it is a real reading of the same equation rather than a failure of it. **Absent in the two cases where there is genuinely no vapour fraction**: * the iteration converged to `x = y = z` (`phase: trivial`), where every `K_i` is 1, `g(beta)` is identically zero, and `beta` is *indeterminate* rather than out of range; and * no Rachford-Rice root exists at all (`phase: all_liquid` or `all_vapour` with `iterations: 1`), which is a proof that the feed is single phase. A caller must branch on `phase` and on presence, never on whether the number looks plausible. The two absences are distinguished by `iterations`: a feed diagnosed without a root settles on the first step, and a trivial solution is the result of an iteration that ran. |
| `x` | dimensionless | liquid-phase mole fractions |
| `y` | dimensionless | vapour-phase mole fractions |
| `k` | dimensionless | K-values, `K_i = y_i / x_i = phi_i^L / phi_i^V`. Reported because they are the iterate the loop converges on, so a caller can see the convergence rather than inferring it from `residual`. |
| `ln_phi_liquid` | dimensionless | fugacity coefficients in the liquid phase, as logarithms. Carried so the identity the iteration drives - that these differ from the vapour ones by `ln K_i` - is checkable from the result alone rather than only from a test. |
| `ln_phi_vapour` | dimensionless | fugacity coefficients in the vapour phase, as logarithms |
| `z_liquid` | dimensionless | the liquid root of the cubic at the liquid composition. The *smallest* admissible root, selected by ordering and never by an initial guess. |
| `z_vapour` | dimensionless | the vapour root at the vapour composition, the largest admissible one |
| `min_t_over_tc` | dimensionless | the smallest `T / Tc_i` over the components - how close the mixture is to the nearest component's critical point. Reported because a bound is checked against it, and a bound whose subject a caller cannot see is a bound that reads as validation while asserting nothing. |
| `phase` | two_phase / all_liquid / all_vapour / trivial | What the converged state is, and the field a caller must branch on. * `two_phase` - `beta` in `[0, 1]`. A genuine split. * `all_liquid` - `beta < 0`, or no Rachford-Rice root exists and every `K_i` is below one. The feed is subcooled liquid. * `all_vapour` - `beta > 1`, or no root exists and every `K_i` is above one. The feed is superheated vapour. * `trivial` - the iteration converged to `x = y = z`. The feed is single phase and **this does not say which one**: the K-values straddled one throughout, so the no-root proof never applied, and distinguishing the two phases is a stability analysis this model does not perform. It is the one outcome that is a statement about the *model* rather than about the feed. All four are reachable and each is covered by a test. The first three are diagnoses of the feed; `trivial` is an admission. |
| `iterations` | dimensionless | successive-substitution steps taken, including the final evaluation of the converged state. Carried because the answer alone does not say whether the loop did any work, and because the cross-language agreement test compares iteration counts as the sharpest cheap check that both implementations ran the same procedure. |
| `residual` | dimensionless | `rms_i |ln K_i - ln K_i_previous|` at the last step the loop completed. Not recomputed after the final evaluation, so it is the quantity the stopping rule actually tested. **NaN when no step completed**, which happens only for a feed diagnosed as single phase without a root - there is nothing to report a residual of, and zero would read as convergence. It is also above the tolerance when the loop stopped on the trivial-solution test rather than on the residual, which is expected rather than a failure: `iterations` and `phase` say which stopping rule fired. |

| Bound | On violation | Why |
|---|---|---|
| `T > 0` | raises | an absolute temperature; zero and below are not states |
| `P > 0` | raises | an absolute pressure; zero and below are not states. This also keeps every `B_i` positive, which `eos.pr_z_factor` requires - `B` is proportional to pressure and the cubic degenerates at zero. |
| `min_t_over_tc <= 0.9` | warns `OUT_OF_VALID_RANGE` | Above `T / Tc = 0.9` for the nearest component the cubic's roots are close to coalescing, the two phases stop being distinguishable, and successive substitution converges slowly or to the trivial solution. This is a warning rather than an error because the arithmetic is still defined and the answer is still the equation's - it is the accuracy and the iteration that degrade, not the meaning. Recorded as a bound rather than left in prose so a caller sees the caveat without having to read this file. |

## Assumptions

- **`Tc`, `Pc`, `omega` and `kij` are the caller's and their correctness is NOT CHECKED.** A databank ships with the library (`azoth.eos.component`), and a caller who supplies their own values replaces it silently - nothing checks that the two agree. A `kij` copied from a table whose convention differs by a sign produces a flash that converges cleanly to a wrong answer; see the same warning in `eos.vdw1f_mix_binary`.
- the equation of state is Peng-Robinson with the coefficient `eos.pr_kappa` computes. PRSV would give different K-values from the same inputs and this model does not accept a coefficient.
- **the mixture fugacity coefficient is not a registered calculation.** `eos.pr_departure` covers a *pure* component; the mixture form, which carries the sum over `x_j a_ij` and the `b_i / b_mix` term, lives in the model layer. That is the one piece of arithmetic here that no kernel checks, and it is mitigated by composition rather than by assertion: at `N = 1` it must reproduce `eos.pr_departure` exactly (the bracket factor collapses to 1), and at `N = 2` the mixture parameters must reproduce `eos.vdw1f_mix_binary` and the vapour fraction `eos.rachford_rice_binary`. Both are tested, in both languages.
- **there is no stability test.** Successive substitution finds a stationary point of the flash equations; whether the feed was stable is a different question, and a converged `trivial` is where the difference shows. Which single phase a `trivial` feed is cannot be answered by this model. This is the largest of the model layer's three: `all_liquid` and `all_vapour` are readings of `beta`, not diagnoses.
- **no damping and no acceleration.** Plain successive substitution, exactly as the scheme names it. Michelsen's minimum-variable Newton is a different scheme and a later milestone; this one converges slowly near the critical point and, as the verification notes record, sometimes to the wrong stationary point.
- the components' `Tc`, `Pc` and `omega` are taken to be mutually consistent and to describe the same substances as the `kij` matrix's indices. None of that is checkable from the numbers.
- the two phases are assumed to be at the same temperature and pressure as the feed - a flash, not a rigorous column stage. Nothing here does an energy balance, and the feed's enthalpy is not an input.

## Cases

| Case | Inputs | Expected |
|---|---|---|
| `methane_and_butane_at_330_k_and_25_bar` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], T = 330.0, P = 2500000.0, z = [0.6, 0.4] | beta = 0.8447220271668119, x = [0.08646861506372552, 0.9135313849362734], y = [0.6943980503344653, 0.3056019496655348], k = [8.030636894354197, 0.3345281341230035], ln_phi_liquid = [2.0769347319092257, -1.4953745844048578], ln_phi_vapour = [-0.006329107265372294, -0.4003402894548284], z_liquid = 0.0933521915309248, z_vapour = 0.8722188526935039, min_t_over_tc = 0.7762514113662025, iterations = 13 |
| `methane_propane_and_butane_at_350_k_and_50_bar` | Tc = [190.56, 369.83, 425.12], Pc = [4599200.0, 4248000.0, 3796000.0], omega = [0.01142, 0.1523, 0.2002], kij = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], T = 350.0, P = 5000000.0, z = [0.5, 0.3, 0.2] | beta = 0.9897493300887433, x = [0.17105357654811876, 0.3837611570605945, 0.445185266391286], y = [0.5034068436348336, 0.2991324995670004, 0.19746065679816602], k = [2.9429775968070517, 0.7794757079069585, 0.44354715150121843], ln_phi_liquid = [1.140973487143631, -0.7934020172498539, -1.631208601116444], ln_phi_vapour = [0.06155163030297489, -0.5442682626142179, -0.8182574350707497], z_liquid = 0.1923272374793064, z_vapour = 0.6962296633271022, min_t_over_tc = 0.8232969514490026, iterations = 22 |
| `a_subcooled_liquid_feed` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], T = 300.0, P = 3000000.0, z = [0.1, 0.9] | beta = -0.046417400850318154, x = [0.1343066534039841, 0.8656933465960167], y = [0.8733970112336372, 0.12660298876638099], k = [6.503006285224948, 0.1462446133659167], ln_phi_liquid = [1.8168271437554724, -2.413916106531643], ln_phi_vapour = [-0.05543743164054027, -0.49144148129974313], z_liquid = 0.1109805757787235, z_vapour = 0.8903679940210696, min_t_over_tc = 0.7056831012420023, iterations = 12 |

## References

- Rachford, H. H.; Rice, J. D. (1952). "Procedure for Use of Electronic Digital Computers in Calculating Flash Vaporization Hydrocarbon Equilibrium." Journal of Petroleum Technology 4(10), 19-3. DOI 10.2118/952327-G. (the vapour-fraction equation; the citation is confirmed and the equation number is not - see eos.rachford_rice_binary, which records how the first draft of that citation was wrong in every part)
- Michelsen, M. L. (1982). "The isothermal flash problem. Part I. Stability." Fluid Phase Equilibria 9(1), 1-19. DOI 10.1016/0378-3812(82)85001-2. (successive substitution with a fugacity-coefficient update, and the stability problem this model does NOT solve; the DOI is confirmed and the paper has not been read)
- Wilson, G. M. (1969). "A Modified Redlich-Kwong Equation of State, Application to General Physical Data Calculations." Paper 15C, AIChE 65th National Meeting. (the K-value estimate the initialisation is named after. The constant 5.373 is quoted as 5.37 in some sources, and the paper is dated 1968 in some and 1969 in others; neither discrepancy is resolved here, and the value this model uses is the one in `initialisation: wilson`.)
