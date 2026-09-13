# Dew-point pressure

`eos.dew_pressure`

The pressure at which a vapour of known composition first condenses, at a fixed temperature. Solved by successive substitution: the pressure, the incipient liquid's composition and the K-values are refined together until the liquid the vapour would produce is exactly one mole.
The companion of `eos.bubble_pressure` and the same iteration with the other phase held. Its equation is `sum_i y_i / K_i = 1` where `y` is the vapour composition that was supplied, and the bubble point's is `sum_i x_i K_i = 1`. The two differ by which of the two compositions is the input, which is why they are two models rather than one with a switch - the argument is called `y` here and `x` there, and a shared name would have to be `z`, which means the feed and would be wrong.

## Source

**Michelsen, M. L. (1982)** (The isothermal flash problem. Part II. Phase-split calculation. Fluid Phase Equilibria 9(1), 21-40. The dew-point iteration is the same successive- substitution scheme as the bubble point's with the roles of the two phases exchanged.
)

DOI: [10.1016/0378-3812(82)85002-4](https://doi.org/10.1016/0378-3812(82)85002-4)

## Notes

`sum_i y_i / K_i = 1` says that the liquid in equilibrium with the vapour is one mole of liquid, and `K_i = phi_i^L / phi_i^V` is the definition of a K-value. What is ours is the initialisation, the update, the tolerance and the guard.

#### The update is the bubble point's, inverted, and getting that wrong is easy

The pressure update is `P <- P / S` here, against the bubble point's `P <- P * S`, and the reason is worth stating because an earlier draft of this model had it the wrong way round and diverged without saying so.

The rule in both models is the same sentence: **if there is too much of the incipient phase, move away from it.** For a bubble point the incipient phase is vapour, `S = sum_i x_i K_i` is how much vapour the liquid would give off, and `S > 1` means too much - so the pressure rises. For a dew point the incipient phase is liquid, `S = sum_i y_i / K_i` is how much liquid the vapour would condense, and `S < 1` means too *little* - so here too the pressure rises, and `P / S` is what raises it.

Measured with the sign wrong: the dew point of methane/n-butane at 300 K ran down to `P = 3.6e-07` Pa over nine iterations and then failed inside `eos.pr_z_factor`, because the pressure had gone to zero rather than to the answer. Nothing about the trajectory was obviously wrong until it was.

#### The trivial solution, and why the guard is on the K-values

`sum_i y_i / K_i = 1` is satisfied by `K_i = 1` for every `i` at *every* pressure, since `sum_i y_i = 1` identically. That trivial solution is a fixed point of the iteration and the physical one is not, so a vapour with no dew point at this temperature converges to it.

Testing the residual does not detect it, and the reason is measured rather than argued: `S - 1 = sum_i y_i (1/K_i - 1)` is a **weighted sum whose terms cancel**. On a state with no dew point the individual deviations can be of order `1e-7` while their weighted sum falls below `1e-12`. So the guard is on the K-values, at `max_i |ln K_i| < 1e-02`. The measurement behind the constant is of the lowest value each state's `max_i |ln K_i|` reaches at any step, not of where it ends up: genuine dew points keep it at `1.30` or above throughout, and states with no dew point drive it below `1e-06`. `eos.bubble_pressure`'s notes carry the full argument and record the drafting error that produced a wrong threshold the first time; it is the same trap and the same constant.

A vapour above the mixture's critical condition has no dew point, and this is how the model says so: `OutOfRange` naming `min_t_over_tc`, the same shape `eos.pure_saturation` uses for `T >= Tc`.

#### No accuracy claim

As everywhere in this namespace, none of the bounds below is an accuracy claim - they are about where the calculation is defined.

## Algorithm

A model rather than a calculation: what this page pins down is the procedure,
not an equation, and both implementations read it from here.

| Setting | Value |
|---|---|
| Scheme | `dew_pressure_successive_substitution` |
| Convergence | `absolute` |
| Tolerance | `1e-12` |
| Max iterations | `200` |
| Initialisation | `wilson_raoult` |

## Inputs

| Name | Unit | Description |
|---|---|---|
| `Tc` | K | critical temperatures, in the mixture's component order |
| `Pc` | Pa | critical pressures, in the same order |
| `omega` | dimensionless | acentric factors, in the same order. All three vectors are the caller's - this library ships a databank of them (`azoth.eos.component`) - and they must be mutually consistent, which is not checked. |
| `kij` | dimensionless | binary interaction parameters, as in `eos.pt_flash` |
| `T` | K | absolute temperature. The dew point is sought at this temperature; the pressure is what is solved for. |
| `y` | dimensionless | the vapour's mole fractions. Checked rather than renormalised, as the flash's feed is. |


## Outputs

| Name | Unit | Description |
|---|---|---|
| `pressure` | Pa | the dew-point pressure |
| `incipient` | dimensionless | the composition of the liquid that first appears - the `x` in equilibrium with the supplied `y`. |
| `k` | dimensionless | K-values at the converged pressure, as in `eos.bubble_pressure` |
| `z_liquid` | dimensionless | the liquid root of the cubic at the converged state |
| `z_vapour` | dimensionless | the vapour root |
| `min_t_over_tc` | dimensionless | the smallest `T / Tc_i` over the components |
| `iterations` | dimensionless | pressure updates taken |
| `residual` | dimensionless | `|sum_i y_i / K_i - 1|` at the last completed step. Not a value a case pins, for the same cancellation reason the guard above exists. |

| Bound | On violation | Why |
|---|---|---|
| `T > 0` | raises | an absolute temperature; zero and below are not states |
| `min_t_over_tc < 1` | warns `OUT_OF_VALID_RANGE` | Above the nearest component's critical temperature the iteration converges more slowly and the answer is less trustworthy. A warning rather than an error because a mixture's critical temperature is not any one component's, so a mixture can still have a dew point where this bound is crossed. |

## Assumptions

- **`Tc`, `Pc`, `omega` and `kij` are the caller's and their correctness is NOT CHECKED.** A databank ships with the library (`azoth.eos.component`), and a caller who supplies their own values replaces it silently - nothing checks that the two agree.
- the equation of state is Peng-Robinson with the coefficient `eos.pr_kappa` computes, and the mixture fugacity coefficient is the one `eos.pt_flash` uses.
- **the model is for mixtures, and refuses a single component.** For one component the dew point is the saturation pressure that `eos.pure_saturation` computes.
- **the pressure is solved for at a fixed temperature.** A dew point at a fixed pressure, with the temperature as the unknown, is a later milestone.
- the vapour composition is the caller's and is taken as given, and its stability is not asked. A vapour that would have condensed earlier is not detected.
- no energy balance. The dew point is a statement about phase equilibrium at one temperature and pressure.

## Cases

| Case | Inputs | Expected |
|---|---|---|
| `methane_and_butane_at_300_k` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], T = 300.0, y = [0.8, 0.2] | pressure = 1568262.871743118, incipient = [0.06533180227306302, 0.934668197726937], k = [12.245184920151232, 0.21397967801460738], z_liquid = 0.059502183160440175, z_vapour = 0.9238233886624134, min_t_over_tc = 0.7056831012420023, iterations = 25 |
| `methane_propane_and_butane_at_300_k` | Tc = [190.56, 369.83, 425.12], Pc = [4599200.0, 4248000.0, 3796000.0], omega = [0.01142, 0.1523, 0.2002], kij = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], T = 300.0, y = [0.5, 0.3, 0.2] | pressure = 1032568.2312776279, incipient = [0.03357294952023628, 0.31034212171707015, 0.6560849287626935], k = [14.892942298646231, 0.9666750950218941, 0.3048385829823517], z_liquid = 0.03812011709239035, z_vapour = 0.902967011675539, min_t_over_tc = 0.7056831012420023, iterations = 20 |

## References

- Michelsen, M. L. (1982). "The isothermal flash problem. Part II. Phase-split calculation." Fluid Phase Equilibria 9(1), 21-40. DOI 10.1016/0378-3812(82)85002-4. (the successive-substitution scheme; the DOI is confirmed and the paper has not been read)
- Wilson, G. M. (1969). "A Modified Redlich-Kwong Equation of State, Application to General Physical Data Calculations." Paper 15C, AIChE 65th National Meeting. (the K-value estimate the initial guess is built from)

## Notes

`sum_i y_i / K_i = 1` says that the liquid in equilibrium with the vapour is one mole of liquid, and `K_i = phi_i^L / phi_i^V` is the definition of a K-value. What is ours is the initialisation, the update, the tolerance and the guard.
# The update is the bubble point's, inverted, and getting that wrong is easy
The pressure update is `P <- P / S` here, against the bubble point's `P <- P * S`, and the reason is worth stating because an earlier draft of this model had it the wrong way round and diverged without saying so.
The rule in both models is the same sentence: **if there is too much of the incipient phase, move away from it.** For a bubble point the incipient phase is vapour, `S = sum_i x_i K_i` is how much vapour the liquid would give off, and `S > 1` means too much - so the pressure rises. For a dew point the incipient phase is liquid, `S = sum_i y_i / K_i` is how much liquid the vapour would condense, and `S < 1` means too *little* - so here too the pressure rises, and `P / S` is what raises it.
Measured with the sign wrong: the dew point of methane/n-butane at 300 K ran down to `P = 3.6e-07` Pa over nine iterations and then failed inside `eos.pr_z_factor`, because the pressure had gone to zero rather than to the answer. Nothing about the trajectory was obviously wrong until it was.
# The trivial solution, and why the guard is on the K-values
`sum_i y_i / K_i = 1` is satisfied by `K_i = 1` for every `i` at *every* pressure, since `sum_i y_i = 1` identically. That trivial solution is a fixed point of the iteration and the physical one is not, so a vapour with no dew point at this temperature converges to it.
Testing the residual does not detect it, and the reason is measured rather than argued: `S - 1 = sum_i y_i (1/K_i - 1)` is a **weighted sum whose terms cancel**. On a state with no dew point the individual deviations can be of order `1e-7` while their weighted sum falls below `1e-12`. So the guard is on the K-values, at `max_i |ln K_i| < 1e-02`. The measurement behind the constant is of the lowest value each state's `max_i |ln K_i|` reaches at any step, not of where it ends up: genuine dew points keep it at `1.30` or above throughout, and states with no dew point drive it below `1e-06`. `eos.bubble_pressure`'s notes carry the full argument and record the drafting error that produced a wrong threshold the first time; it is the same trap and the same constant.
A vapour above the mixture's critical condition has no dew point, and this is how the model says so: `OutOfRange` naming `min_t_over_tc`, the same shape `eos.pure_saturation` uses for `T >= Tc`.
# No accuracy claim
As everywhere in this namespace, none of the bounds below is an accuracy claim - they are about where the calculation is defined.
