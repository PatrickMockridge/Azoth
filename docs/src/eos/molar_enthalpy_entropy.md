# Molar enthalpy and entropy of a mixture

`eos.molar_enthalpy_entropy`

The absolute molar enthalpy and entropy of a mixture at a state, composed from three parts: the ideal-gas enthalpy and entropy carried in from a reference state, the ideal-gas heat capacity's integral between that reference and the state, and the departure functions for the equation of state.
``` H = H_ig(T_ref) + integral Cp dT      + H_dep S = S_ig(T_ref) + integral Cp/T dT    - R ln(P/P_ref) - R sum z_i ln z_i + S_dep ```
A *direct* model rather than a procedure: there is no iteration here, and no algorithm to pin down. It is a model rather than a calculation because its inputs are composition vectors, which the calc registry - whose inputs are scalars - has nowhere to put, and because a library function with one implementation would be the one thing this project is organised against. See `kind` below.

## Source

**The composition of ideal-gas and residual contributions** (Standard thermodynamics: an absolute enthalpy or entropy is a reference value plus an ideal-gas integral plus a departure, and the split between them is a convention rather than a result. No source states this particular assembly, because the assembly is not a discovery - what matters is that the terms are written down and the datum is the caller's, which is what this spec does.
)

## Notes

Nothing here is a correlation. The departure functions are `eos.pr_departure`'s machinery applied to a mixture, the ideal-gas integrals are `eos.ideal_gas_cp`'s polynomial integrated exactly, and the mixing term is the ideal-gas entropy of mixing. What is ours, and what this entry is for, is the **assembly** - which terms belong where, and what the caller is responsible for.

#### The datum is the caller's, and that is the dangerous part

An absolute enthalpy is absolute only relative to a reference state, and this library ships no reference. `h_ref` and `s_ref` are the caller's - one per component, the ideal-gas values at `T_ref` and `P_ref` - and nothing checks that they come from the same source, or that two calls using different ones were meant to be compared.

**Two enthalpies from different datums are not comparable, and subtracting them produces a plausible number rather than an error.** A heat of reaction computed from two streams whose `h_ref` came from different tables is wrong by the difference between the tables' reference states, which for a formation-property datum can be hundreds of kJ/mol. This is the next silent wrong number this programme would otherwise produce, and there is no arithmetic that catches it - so it is stated here, in the API docstring, and in `assumptions`, and the worked example uses a zero datum so that a reader comparing it against their own numbers knows exactly which datum it is on.

#### Why the ideal-gas entropy has two terms that the enthalpy does not

With every coefficient and reference value set to zero, `H` is **exactly** the departure - `R*T*h_dep_rt`, which is what `eos.pr_departure` computes for a pure component and what the mixture module computes for a mixture. `S` is not, and the difference is structural rather than an oversight: an ideal gas's entropy depends on pressure and on mixing even when its heat capacity is zero, so

```text S(coefficients = 0, refs = 0) = -R ln(P/P_ref) - R sum_i z_i ln z_i + R*s_dep_r ```

and the two leading terms are ideal-gas terms that no coefficient switches off. The spec case below pins that, because a reader who expected `S = S_dep` exactly would take the difference for a defect.

#### The mixing term is inside the ideal-gas part, not the departure

`-R sum_i z_i ln z_i` is added explicitly and is **not** part of `s_dep_r`. The departure functions are defined against the ideal-gas *mixture* at the same temperature, pressure and composition - which is why `sum_i z_i ln phi_i` vanishes for an ideal mixture - so the entropy of mixing belongs to the ideal-gas side. A reader who expected it in the departure would find the two disagreeing by exactly that term.

#### The composition is a pure function of the state, and the root is the caller's

The model takes the compressibility factor as an input rather than solving for it. That is deliberate: which root describes the phase is a *choice*, a caller holding a `Z` from `eos.pr_z_factor` or from a flash has already made it, and a model that re-derived it would silently overrule them. It also keeps this model free of the root-selection question that `eos.pt_flash` had to answer with a phase enum.

#### No accuracy claim

The departure functions inherit Peng-Robinson's accuracy; the ideal-gas integrals inherit the coefficients'. No bound below is an accuracy claim.

## What this model is

A **direct** model: a computation over vectors, with no iteration and therefore no algorithm block. It is here because a calculation's inputs are scalars, so a composition vector has nowhere in the calc registry to go - and it is a model rather than a library function so that it has two implementations like everything else.

## Inputs

| Name | Unit | Description |
|---|---|---|
| `Tc` | K | critical temperatures, in the mixture's component order |
| `Pc` | Pa | critical pressures, in the same order |
| `omega` | dimensionless | acentric factors, in the same order |
| `kij` | dimensionless | binary interaction parameters, as in `eos.pt_flash` |
| `cp_a` | dimensionless | the constant term of each component's `Cp/R` polynomial. As in `eos.ideal_gas_cp`, the four coefficients are dimensionless because the polynomial is divided through by `R` and written against `T/(1000 K)`. |
| `cp_b` | dimensionless | the coefficient of `theta` in each component's polynomial |
| `cp_c` | dimensionless | the coefficient of `theta**2` in each component's polynomial |
| `cp_d` | dimensionless | the coefficient of `theta**3` in each component's polynomial |
| `h_ref` | J/mol | each component's ideal-gas molar enthalpy at `T_ref`, which is the caller's datum. **Nothing checks it against any other source**, and two enthalpies computed from different datums are not comparable. |
| `s_ref` | J/(mol*K) | each component's ideal-gas molar entropy at `T_ref` and `P_ref` |
| `T_ref` | K | the temperature the reference values are given at |
| `P_ref` | Pa | the pressure `s_ref` is given at. It does not enter the enthalpy. |
| `T` | K | absolute temperature of the state wanted |
| `P` | Pa | absolute pressure of the state wanted |
| `z` | dimensionless | the mixture's mole fractions, checked rather than renormalised |
| `compressibility` | dimensionless | the cubic's root for this phase, `Z = P v /(R T)`. An input rather than something this model solves for: which root describes the phase is the caller's choice, and a caller holding a `Z` from `eos.pr_z_factor` or from a flash has already made it. |


## Outputs

| Name | Unit | Description |
|---|---|---|
| `h` | J/mol | the molar enthalpy, `h_ideal + h_departure` |
| `s` | J/(mol*K) | the molar entropy, `s_ideal + s_departure` |
| `h_ideal` | J/mol | the ideal-gas part: the reference enthalpies plus the heat-capacity integrals. Reported separately because it is the term that carries the datum, so a caller can see which part of the answer came from where. |
| `s_ideal` | J/(mol*K) | the ideal-gas part of the entropy - the reference values, the heat-capacity integrals, `-R ln(P/P_ref)` for the pressure and `-R sum_i z_i ln z_i` for the mixing. The last two do not vanish when the heat-capacity coefficients do. |
| `h_departure` | J/mol | the residual enthalpy, `R*T*h_dep_rt` |
| `s_departure` | J/(mol*K) | the residual entropy, `R*s_dep_r`. **Does not include the entropy of mixing**, which is an ideal-gas term; see the notes. |
| `psi_bar` | dimensionless | the composition-weighted average of the components' `psi` that the mixture departure functions are built on. Reported because it is what makes the mixture departure differ from a pure component's, and because at `N = 1` it must equal that component's own `psi` exactly - a test asserts it. |

| Bound | On violation | Why |
|---|---|---|
| `T > 0` | raises | an absolute temperature; zero and below are not states |
| `T_ref > 0` | raises | the reference temperature, and the argument of `ln(T/T_ref)` in the entropy integral - so zero is a division rather than a state |
| `P > 0` | raises | an absolute pressure; zero and below are not states |
| `P_ref > 0` | raises | the reference pressure, and the argument of `ln(P/P_ref)` |
| `compressibility > 0` | raises | `Z = P v /(R T)` is positive for any state with a positive volume. The stricter check that actually matters - that the root exceeds the mixture's `B`, without which `ln(Z - B)` is the logarithm of a negative number - needs the composition and belongs to the implementation rather than to this table. |

## Assumptions

- **`h_ref` and `s_ref` are the caller's datum, and their correctness and consistency are NOT CHECKED.** Two enthalpies computed from different datums are not comparable, and subtracting them gives a plausible number rather than an error. This is the single most likely way to get a wrong answer out of this model.
- **`Tc`, `Pc`, `omega`, `kij` and the four heat-capacity coefficients are the caller's, and their correctness is NOT CHECKED.** A databank of critical constants ships with the library; the heat-capacity coefficients are not in it.
- the heat-capacity polynomial is valid over the range it was fitted to, and `T`, `T_ref` and the integration path between them all lie in it. The coefficients describe the ideal-gas heat capacity at every temperature on that path, which is an assumption the caller makes and this model cannot.
- the departure functions are relative to the ideal-gas *mixture* at the same temperature, pressure and composition, not to the pure ideal gases. The entropy of mixing is therefore in `s_ideal` and not in `s_departure`.
- **which root `compressibility` names is the caller's choice and is not checked against the phase it describes.** Passing a liquid root for a vapour state gives a number that is internally consistent and describes a state that does not exist.
- no energy balance and no reaction. This is a property of one state, and what a caller does with two of them - a difference, a heat of mixing, a reaction enthalpy - is a composition this model does not make and cannot check.

## Cases

| Case | Inputs | Expected |
|---|---|---|
| `methane_and_butane_with_a_zero_datum` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], cp_a = [4.0, 4.0], cp_b = [1.0, 1.0], cp_c = [-0.5, -0.5], cp_d = [0.1, 0.1], h_ref = [0.0, 0.0], s_ref = [0.0, 0.0], T_ref = 298.15, P_ref = 101325.0, T = 330.0, P = 2500000.0, z = [0.6, 0.4], compressibility = 0.8274482588400789 | h = -351.24593942874526, s = -20.553517170090107, h_ideal = 1130.1847399375104, s_ideal = -17.456668341283496, h_departure = -1481.4306793662556, s_departure = -3.0968488288066105, psi_bar = -0.5619363510353611 |
| `the_same_state_with_the_ideal_gas_terms_off` | Tc = [190.56, 425.12], Pc = [4599200.0, 3796000.0], omega = [0.01142, 0.2002], kij = [[0.0, 0.05], [0.05, 0.0]], cp_a = [0.0, 0.0], cp_b = [0.0, 0.0], cp_c = [0.0, 0.0], cp_d = [0.0, 0.0], h_ref = [0.0, 0.0], s_ref = [0.0, 0.0], T_ref = 298.15, P_ref = 101325.0, T = 330.0, P = 2500000.0, z = [0.6, 0.4], compressibility = 0.8274482588400789 | h = -1481.4306793662556, s = -24.154898040804962, h_ideal = 0.0, s_ideal = -21.058049211998352, h_departure = -1481.4306793662556, s_departure = -3.0968488288066105, psi_bar = -0.5619363510353611 |

## References

- Poling, B. E.; Prausnitz, J. M.; O'Connell, J. P. "The Properties of Gases and Liquids", 5th ed. (the split of an absolute enthalpy into a reference value, an ideal-gas integral and a departure; the edition is described rather than page-cited, and is unconfirmed)
- Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State." Industrial & Engineering Chemistry Fundamentals 15(1), 59-64. DOI 10.1021/i160057a011. (the equation of state the departure functions come from; the citation is confirmed and the equation numbers are not)

## Notes

Nothing here is a correlation. The departure functions are `eos.pr_departure`'s machinery applied to a mixture, the ideal-gas integrals are `eos.ideal_gas_cp`'s polynomial integrated exactly, and the mixing term is the ideal-gas entropy of mixing. What is ours, and what this entry is for, is the **assembly** - which terms belong where, and what the caller is responsible for.
# The datum is the caller's, and that is the dangerous part
An absolute enthalpy is absolute only relative to a reference state, and this library ships no reference. `h_ref` and `s_ref` are the caller's - one per component, the ideal-gas values at `T_ref` and `P_ref` - and nothing checks that they come from the same source, or that two calls using different ones were meant to be compared.
**Two enthalpies from different datums are not comparable, and subtracting them produces a plausible number rather than an error.** A heat of reaction computed from two streams whose `h_ref` came from different tables is wrong by the difference between the tables' reference states, which for a formation-property datum can be hundreds of kJ/mol. This is the next silent wrong number this programme would otherwise produce, and there is no arithmetic that catches it - so it is stated here, in the API docstring, and in `assumptions`, and the worked example uses a zero datum so that a reader comparing it against their own numbers knows exactly which datum it is on.
# Why the ideal-gas entropy has two terms that the enthalpy does not
With every coefficient and reference value set to zero, `H` is **exactly** the departure - `R*T*h_dep_rt`, which is what `eos.pr_departure` computes for a pure component and what the mixture module computes for a mixture. `S` is not, and the difference is structural rather than an oversight: an ideal gas's entropy depends on pressure and on mixing even when its heat capacity is zero, so
```text S(coefficients = 0, refs = 0) = -R ln(P/P_ref) - R sum_i z_i ln z_i + R*s_dep_r ```
and the two leading terms are ideal-gas terms that no coefficient switches off. The spec case below pins that, because a reader who expected `S = S_dep` exactly would take the difference for a defect.
# The mixing term is inside the ideal-gas part, not the departure
`-R sum_i z_i ln z_i` is added explicitly and is **not** part of `s_dep_r`. The departure functions are defined against the ideal-gas *mixture* at the same temperature, pressure and composition - which is why `sum_i z_i ln phi_i` vanishes for an ideal mixture - so the entropy of mixing belongs to the ideal-gas side. A reader who expected it in the departure would find the two disagreeing by exactly that term.
# The composition is a pure function of the state, and the root is the caller's
The model takes the compressibility factor as an input rather than solving for it. That is deliberate: which root describes the phase is a *choice*, a caller holding a `Z` from `eos.pr_z_factor` or from a flash has already made it, and a model that re-derived it would silently overrule them. It also keeps this model free of the root-selection question that `eos.pt_flash` had to answer with a phase enum.
# No accuracy claim
The departure functions inherit Peng-Robinson's accuracy; the ideal-gas integrals inherit the coefficients'. No bound below is an accuracy claim.
