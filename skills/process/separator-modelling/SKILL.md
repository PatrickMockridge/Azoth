# Separator modelling

Educational separator screening indicators for gas load, residence time, and capacity warnings. USE WHEN: a task needs public separator capacity screening without proprietary design methods.

This is a `screening` skill and **partly ported**. azoth's `process.separator` is an
**equilibrium flash split**, and given a vessel's internal diameter and design K it also
reports that vessel's steady gas capacity: the Souders-Brown **utilisation**, the maximum
allowable and actual superficial **gas velocities**, the **gas load factor** and an
**overload flag**. Residence times, droplet cut size, inlet momentum, the droplet-performance
calculator, mesh-pad sizing and the mechanical design are still not ported.

## When to Use

Educational separator screening indicators for gas load, residence time, and capacity warnings. USE WHEN: a task needs public separator capacity screening without proprietary design methods.

## Inputs

The split takes the feed's record (`components`, `feed_n`, `feed_z`, `feed_p`, `feed_t`), a
`pressure_drop`, an entrainment `gas_in_liquid` and an optional `heat_input`. The capacity
surface adds four optional inputs:

- `internal_diameter` (m) and `design_gas_load_factor` (m/s), **both or neither**;
- `orientation` (`horizontal` or `vertical`) and `design_liquid_level_fraction`
  (dimensionless, default `0.8`), which shape the gas area and are refused without the
  diameter and load factor.

## Outputs

The two outlet records (`vapour_*`, `liquid_*`) and, when the vessel is stated:
`capacity_utilization`, `max_allowable_gas_velocity`, `gas_superficial_velocity`,
`gas_load_factor` and `is_overloaded`. The three velocities are absent when the system
carries no gas; a liquid-only state answers utilisation `0.0` and overload `false`.

## How a calculation runs

The feed is flashed at `feed_p - pressure_drop` and at the feed's own temperature (or at the
enthalpy a `heat_input` implies). The capacity surface reads the system that flash leaves:
the gas's volumetric flow over the design gas area, the Souders-Brown velocity from the
design K, and the actual load factor from the operating velocity.

**Two density rules, deliberately separate.** The utilisation and maximum velocity use the
liquid's own density, or `1000.0` where there is no liquid, with no floor. The gas load
factor applies NeqSim's `10.0` kg/m3 near-dry floor first. A near-dry state makes the two
disagree by construction.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
result = azoth.process.separator(
    ["methane", "n-butane"],
    q(1.0, "mol/s"),
    [0.7, 0.3],
    q(20.0, "bar"),
    q(300.0, "K"),
    q(0.0, "Pa"),
    0.0,
    internal_diameter=q(1.0, "m"),
    design_gas_load_factor=q(0.11, "m/s"),
    orientation="vertical",
)
print(result.capacity_utilization, result.gas_load_factor, result.is_overloaded)
```

## The keycard

Not applicable — the fluid is the caller's composition and the databank resolves it.

## Validation checklist

- [ ] Confirm the task actually needs this calculation.
- [ ] State **both** the diameter and the design K, and the orientation the vessel really has.
- [ ] Treat the result as a screening indicator; use NeqSim's validated engine for design.
- [ ] Check whether a single-phase state is being read: a dry gas uses the `1000.0` fallback.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A utilisation five times too low | A vertical scrubber left on the horizontal default (`0.8` derates the area) | State `orientation="vertical"` |
| `gas_load_factor` used as the design K | It is the **actual** Souders-Brown factor | `design_gas_load_factor` is the vessel's basis |
| A horizontal vessel `max_allowable_gas_velocity` changed by `design_liquid_level_fraction` | The velocity is density-based; only the area derates | Compare `max_allowable_gas_flow` instead |
| Half a vessel quietly ignored | Only one of diameter and K stated | Both are required; the port refuses one alone |

## Limitations

One condensed phase only. The class's `SeparatorCapacityAssessment`, `getGasLoadFactor(int)`,
`getDeRatedGasLoadFactor`, the droplet calculator, `sizeFromFlow`/`autoSize` and the
mechanical design are not ported. Residence times are not computed.

## Related Azoth functionality

`process.separator`'s capacity surface (this skill's utilisation, velocities and load
factor), and `process.gas_scrubber`, whose stream side is the same flash and whose own
capacity metric is the utilisation alone.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
